################################################################################
# ORF Analysis - Unified Pipeline Integration
#
# This script ENRICHES blast.parquet with ORF analysis results rather than
# creating separate output files. This enables a single unified pipeline.
#
# Behavior:
#   - Reads blast.parquet
#   - For each species, finds sequences containing ORFs (>= min_orf_length)
#   - Adds columns: ORF_Filtered (bool), ORF_COUNT (int), ORF_SUMMARY (string)
#   - Original rows get ORF_Filtered=FALSE
#   - ORF-passing rows get ORF_Filtered=TRUE with ORF metadata
#   - Writes enriched data back to blast.parquet/blast.csv
#   - Creates orf_analysis.flag marker file
#   - Checks for downstream outputs and warns if they need rerunning
#   - Exports per-species FASTA files (for direct use or debugging)
#
# Usage:
#   Rscript orf_analyser.R <blast_parquet> <species_db> <output_folder> [min_orf_length]
#
# Arguments:
#   blast_parquet  - Path to blast.parquet file (will be modified in place)
#   species_db     - Base path to genome FASTA files (species named {Species}.fa)
#   output_folder  - Output directory for FASTA files and flag
#   min_orf_length - Minimum ORF length in nucleotides (default: 200)
################################################################################

# ==============================================================================
# 0) DEPENDENCIES
# ==============================================================================
options(warn = -1)  # silence warnings
suppressMessages({
  library(tidyverse)     # data wrangling, readr for CSV
  library(arrow)         # read/write parquet
  library(GenomicRanges) # GRanges container and helpers
  library(plyranges)     # dplyr-like verbs for GRanges
  library(rtracklayer)   # GFF3 export
  library(Biostrings)    # FASTA IO, getSeq, reverse-complement handling
  library(ORFik)         # ORF finding on DNAStringSet / GRanges
  library(IRanges)       # IRanges utils: elementNROWS, start/end on ranges
})

# ==============================================================================
# 1) COMMAND-LINE ARGUMENTS
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript orf_analyser.R <blast_parquet> <species_db> <output_folder> [min_orf_length]")
}

args.blast_parquet   <- args[1]  # Path to blast.parquet
args.species_db      <- args[2]  # Base path to genome FASTA files
args.output_folder   <- args[3]  # Output dir for FASTA files and flag
args.min_orf_length  <- if (length(args) >= 4) as.integer(args[4]) else 200L

# Derive paths
blast_dir <- dirname(args.blast_parquet)
blast_csv_path <- file.path(blast_dir, "blast.csv")
flag_file_path <- file.path(args.output_folder, "orf_analysis.flag")

# Create output directories
dir.create(args.output_folder, showWarnings = FALSE, recursive = TRUE)
fasta_output_dir <- file.path(args.output_folder, "fastas")
dir.create(fasta_output_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting ORF analysis (unified pipeline mode)...")
message("  BLAST parquet: ", args.blast_parquet)
message("  Species DB:    ", args.species_db)
message("  Output folder: ", args.output_folder)
message("  Min ORF len:   ", args.min_orf_length)

# ==============================================================================
# 2) CHECK FOR EXISTING DOWNSTREAM OUTPUTS (WARNING SYSTEM)
# ==============================================================================
check_downstream_outputs <- function(blast_dir) {
  downstream_files <- c(
    file.path(blast_dir, "rpsblast.parquet"),
    file.path(blast_dir, "rpsblast.csv"),
    file.path(blast_dir, "domains.parquet"),
    file.path(blast_dir, "domains.csv")
  )

  existing <- downstream_files[file.exists(downstream_files)]

  if (length(existing) > 0) {
    message("\n")
    message("==============================================================")
    message("WARNING: Downstream outputs already exist!")
    message("==============================================================")
    message("The following files were generated before this ORF analysis:")
    for (f in existing) {
      message("  - ", basename(f))
    }
    message("")
    message("These files were generated using ALL BLAST hits.")
    message("After this ORF analysis completes, downstream rules will use")
    message("only ORF-filtered sequences. You should rerun:")
    message("  --rpsblast --rpsbproc --rpsbproc_parser")
    message("==============================================================\n")
  }
}

check_downstream_outputs(blast_dir)

# ==============================================================================
# 3) HELPER FUNCTIONS
# ==============================================================================

# Normalizer for IDs (vectorized)
normalize_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)                         # remove surrounding whitespace
  x <- sub("^[^|]*\\|", "", x)           # drop leading DB prefix like 'ref|'/'gb|'/'emb|'
  x <- sub("\\|$", "", x)                # drop trailing pipe if present
  x <- sub("^chr", "", x, ignore.case = TRUE) # drop 'chr' prefix (case-insensitive)
  x <- sub("\\.\\d+$", "", x)            # drop version suffix '.number' at end
  x
}

# Summarize ORFs per slice for CSV
summarize_orfs <- function(ir, parent_gr) {
  # ir: IRanges of ORFs in local coordinates (relative to the extracted slice)
  # parent_gr: the GRanges range for that slice (genome coordinates + strand)
  if (length(ir) == 0) return(NA_character_)  # safeguard

  ps <- IRanges::start(ir)  # local starts
  pe <- IRanges::end(ir)    # local ends

  if (as.character(strand(parent_gr)) == "+") {
    # '+' strand: shift local to genome by adding (slice_start - 1)
    gs <- start(parent_gr) + ps - 1L
    ge <- start(parent_gr) + pe - 1L
  } else {
    # '-' strand: mirror around slice_end
    gs <- end(parent_gr) - pe + 1L
    ge <- end(parent_gr) - ps + 1L
  }

  # Rough amino-acid length: floor(nt_len / 3)
  aa <- floor((pe - ps + 1L) / 3L)

  # Build items like "seqid:123-456(-):111aa"
  items <- paste0(
    as.character(seqnames(parent_gr)), ":", gs, "-", ge,
    "(", as.character(strand(parent_gr)), "):", aa, "aa"
  )

  # Join multiple ORFs with ';'
  paste(items, collapse = ";")
}

# Process a single species - returns data frame with ORF columns added
process_species_for_enrichment <- function(species, d0, args.species_db,
                                            fasta_output_dir, args.min_orf_length) {

  message(sprintf("Processing species: %s", species))

  # Build genome FASTA path
  genome_fasta <- file.path(args.species_db, paste0(species, ".fa"))

  # Check if genome FASTA exists
  if (!file.exists(genome_fasta)) {
    warning(sprintf("Genome FASTA not found for %s at %s. Skipping.", species, genome_fasta))
    return(list(orf_data = NULL, fasta_exported = FALSE))
  }

  # Filter BLAST hits to this species
  d1 <- d0 %>%
    filter(Species == species) %>%
    filter(!is.na(`Subject ID`), !is.na(`S. Start`), !is.na(`S. End`)) %>%
    mutate(
      `S. Start` = as.integer(`S. Start`),
      `S. End`   = as.integer(`S. End`)
    )

  if (nrow(d1) == 0) {
    warning(sprintf("No BLAST hits for species %s. Skipping.", species))
    return(list(orf_data = NULL, fasta_exported = FALSE))
  }

  message(sprintf("  Loaded %d hits for %s", nrow(d1), species))

  # Strand inference + coordinate normalization
  d2 <- d1 %>%
    mutate(
      strand = if_else(`S. Start` < `S. End`, "+", "-"),
      start  = pmin(`S. Start`, `S. End`),
      end    = pmax(`S. Start`, `S. End`)
    )

  # Build GRanges from the data frame
  gr <- makeGRangesFromDataFrame(
    d2,
    seqnames.field      = "Subject ID",
    start.field         = "start",
    end.field           = "end",
    strand.field        = "strand",
    keep.extra.columns  = TRUE
  )

  # Read genome FASTA
  genome <- tryCatch(
    readDNAStringSet(genome_fasta),
    error = function(e) {
      warning(sprintf("Error reading genome FASTA for %s: %s", species, e$message))
      return(NULL)
    }
  )

  if (is.null(genome)) return(list(orf_data = NULL, fasta_exported = FALSE))

  # Harmonize sequence identifiers
  gn_raw         <- names(genome)
  gn_norm        <- normalize_ids(gn_raw)
  seq_lvls       <- seqlevels(gr)
  seq_lvls_norm  <- normalize_ids(seq_lvls)

  # Map GRanges levels to FASTA names
  mapped <- gn_raw[match(seq_lvls_norm, gn_norm)]
  ok     <- !is.na(mapped)

  # Rename matched levels
  if (any(ok)) {
    gr <- renameSeqlevels(gr, setNames(mapped[ok], seq_lvls[ok]))
  }

  # Keep only seqlevels present in FASTA
  keep <- intersect(seqlevels(gr), names(genome))
  drop <- setdiff(seqlevels(gr), keep)
  if (length(drop)) {
    message(sprintf("  Dropping %d unmapped seqlevels for %s", length(drop), species))
  }
  gr <- keepSeqlevels(gr, keep, pruning.mode = "coarse")

  if (length(gr) == 0) {
    warning(sprintf("No valid ranges after seqlevel filtering for %s. Skipping.", species))
    return(list(orf_data = NULL, fasta_exported = FALSE))
  }

  # Extract sequences for each range (strand-aware)
  seqs <- getSeq(genome, gr)

  # Find ORFs in each extracted sequence
  orfs <- ORFik::findORFs(
    seqs,
    longestORF     = TRUE,
    minimumLength  = args.min_orf_length
  )

  # Build indices of sequences with >= 1 ORF
  has_idx <- if (length(orfs) == length(seqs)) {
    which(IRanges::elementNROWS(orfs) > 0)
  } else {
    sort(unique(as.integer(names(orfs))))
  }

  if (length(has_idx) == 0) {
    message(sprintf("  No sequences with ORFs found for %s.", species))
    return(list(orf_data = NULL, fasta_exported = FALSE))
  }

  message(sprintf("  Found %d sequences with ORFs (of %d total) for %s",
                  length(has_idx), length(gr), species))

  # Filter GRanges and data frame to those with >= 1 ORF
  gr_filt <- gr[has_idx]
  d2_filt <- d2[has_idx, , drop = FALSE]

  # Align ORF list to filtered GRanges
  orfs_filt <- if (length(orfs) == length(seqs)) {
    orfs[has_idx]
  } else {
    orfs[match(has_idx, as.integer(names(orfs)))]
  }

  # Add ORF columns to filtered data frame
  d2_filt$ORF_COUNT <- as.integer(IRanges::elementNROWS(orfs_filt))
  d2_filt$ORF_SUMMARY <- vapply(
    seq_along(gr_filt),
    function(j) summarize_orfs(orfs_filt[[j]], gr_filt[j]),
    FUN.VALUE = character(1)
  )
  d2_filt$ORF_Filtered <- TRUE

  # Remove helper columns (strand, start, end) that we added
  d2_filt <- d2_filt %>%
    select(-any_of(c("strand", "start", "end")))

  # Build tag vector for IDs
  tag_vec <- mcols(gr_filt)$Tag
  if (is.null(tag_vec)) tag_vec <- NA_character_
  tag_vec <- ifelse(is.na(tag_vec) | tag_vec == "", "NA", as.character(tag_vec))

  # === Export FASTA (for backward compatibility and debugging) ===
  fa <- getSeq(genome, gr_filt)

  # FASTA headers: <seqid>:<start>-<end>|tag:<tag>
  # Match the format expected by downstream tools
  hdr <- paste0(
    as.character(seqnames(gr_filt)), ":",
    start(gr_filt), "-", end(gr_filt), "|tag:",
    tag_vec
  )
  names(fa) <- hdr

  fasta_path <- file.path(fasta_output_dir, paste0(species, ".fa"))
  writeXStringSet(fa, fasta_path)
  message(sprintf("  Wrote FASTA: %s", fasta_path))

  return(list(
    orf_data = d2_filt,
    fasta_exported = TRUE,
    total_hits = nrow(d1),
    hits_with_orfs = length(has_idx)
  ))
}

# ==============================================================================
# 4) LOAD BLAST TABLE
# ==============================================================================
message("Loading BLAST data...")
d0 <- arrow::read_parquet(args.blast_parquet)

# Check if ORF analysis was already run
if ("ORF_Filtered" %in% colnames(d0)) {
  message("WARNING: ORF_Filtered column already exists in blast.parquet")
  message("Removing previous ORF analysis results before re-running...")
  d0 <- d0 %>%
    filter(ORF_Filtered == FALSE | is.na(ORF_Filtered)) %>%
    select(-any_of(c("ORF_Filtered", "ORF_COUNT", "ORF_SUMMARY")))
}

species_list <- unique(d0$Species)
message(sprintf("Found %d species in BLAST parquet", length(species_list)))

# ==============================================================================
# 5) PROCESS ALL SPECIES
# ==============================================================================
orf_results <- lapply(species_list, function(species) {
  tryCatch(
    process_species_for_enrichment(species, d0, args.species_db,
                                   fasta_output_dir, args.min_orf_length),
    error = function(e) {
      warning(sprintf("Error processing %s: %s", species, e$message))
      return(list(orf_data = NULL, fasta_exported = FALSE))
    }
  )
})

# ==============================================================================
# 6) COMBINE ORF-FILTERED DATA
# ==============================================================================
# Collect all ORF-filtered data frames
orf_data_list <- lapply(orf_results, function(x) x$orf_data)
orf_data_list <- Filter(Negate(is.null), orf_data_list)

if (length(orf_data_list) == 0) {
  message("\nNo species produced ORF-filtered sequences.")
  message("blast.parquet will not be modified.")
  # Still create the flag file to indicate ORF analysis was attempted
  writeLines(
    c(
      paste0("timestamp: ", Sys.time()),
      paste0("min_orf_length: ", args.min_orf_length),
      "status: no_orfs_found"
    ),
    flag_file_path
  )
  message(sprintf("Wrote flag file: %s", flag_file_path))
  quit(save = "no", status = 0)
}

orf_combined <- bind_rows(orf_data_list)
message(sprintf("\nCombined %d ORF-filtered sequences from %d species",
                nrow(orf_combined), length(orf_data_list)))

# ==============================================================================
# 7) MERGE WITH ORIGINAL DATA
# ==============================================================================
# Add ORF columns to original data (all FALSE)
d0_enriched <- d0 %>%
  mutate(
    ORF_Filtered = FALSE,
    ORF_COUNT = NA_integer_,
    ORF_SUMMARY = NA_character_
  )

# Ensure column order matches
common_cols <- intersect(colnames(d0_enriched), colnames(orf_combined))
d0_enriched <- d0_enriched %>% select(all_of(common_cols))
orf_combined <- orf_combined %>% select(all_of(common_cols))

# Combine: original (ORF_Filtered=FALSE) + ORF-passing (ORF_Filtered=TRUE)
blast_enriched <- bind_rows(d0_enriched, orf_combined)

message(sprintf("Final enriched table: %d rows (%d original + %d ORF-filtered)",
                nrow(blast_enriched), nrow(d0_enriched), nrow(orf_combined)))

# ==============================================================================
# 8) WRITE ENRICHED DATA BACK
# ==============================================================================
message("Writing enriched blast.parquet...")
arrow::write_parquet(blast_enriched, args.blast_parquet)
message(sprintf("  Updated: %s", args.blast_parquet))

message("Writing enriched blast.csv...")
readr::write_csv(blast_enriched, blast_csv_path)
message(sprintf("  Updated: %s", blast_csv_path))

# ==============================================================================
# 9) WRITE FLAG FILE AND MANIFEST
# ==============================================================================
# Flag file indicates ORF analysis has been run
writeLines(
  c(
    paste0("timestamp: ", Sys.time()),
    paste0("min_orf_length: ", args.min_orf_length),
    paste0("species_count: ", length(orf_data_list)),
    paste0("orf_filtered_sequences: ", nrow(orf_combined)),
    "status: complete"
  ),
  flag_file_path
)
message(sprintf("Wrote flag file: %s", flag_file_path))

# Write species manifest (species that have ORF-filtered sequences)
species_with_orfs <- unique(orf_combined$Species)
manifest_path <- file.path(args.output_folder, "species_manifest.txt")
writeLines(species_with_orfs, manifest_path)
message(sprintf("Wrote species manifest: %s (%d species)",
                manifest_path, length(species_with_orfs)))

# ==============================================================================
# 10) SUMMARY
# ==============================================================================
stats <- lapply(orf_results, function(x) {
  if (!is.null(x$orf_data)) {
    list(total = x$total_hits, with_orfs = x$hits_with_orfs)
  } else {
    NULL
  }
})
stats <- Filter(Negate(is.null), stats)

message("\n=== ORF Analysis Summary ===")
message(sprintf("Species with ORFs: %d of %d", length(stats), length(species_list)))
message(sprintf("Total hits analyzed: %d", sum(sapply(stats, function(x) x$total))))
message(sprintf("Total hits with ORFs: %d", sum(sapply(stats, function(x) x$with_orfs))))
message("")
message("Pipeline behavior after ORF analysis:")
message("  - blast_parser will export only ORF_Filtered=TRUE sequences to FASTA")
message("  - Downstream rules (rpsblast, rpsbproc) will process ORF-filtered data")
message("")
message("Done.")
