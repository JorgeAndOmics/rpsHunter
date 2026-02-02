################################################################################
# ORF-aware extraction and export from BLAST hits
# - Loads BLAST hits (parquet)
# - Infers strand and normalizes coordinates
# - Builds GRanges and harmonizes sequence IDs to FASTA headers
# - Extracts sequences honoring strand
# - Calls ORFs and keeps only slices containing >= 1 ORF
# - Exports per-species: CSV (with ORF_COUNT and ORF_SUMMARY), GFF3, and FASTA
#
# Usage:
#   Rscript orf_analyser.R <blast_parquet> <species_db> <output_folder> <min_orf_length>
#
# Arguments:
#   blast_parquet  - Path to blast.parquet file
#   species_db     - Base path to genome FASTA files (species named {Species}.fa)
#   output_folder  - Output directory for all ORF outputs (fastas/, *.csv, *.gff3)
#   min_orf_length - Minimum ORF length in nucleotides (default: 200)
################################################################################

# ==============================================================================
# 0) DEPENDENCIES
# ==============================================================================
options(warn = -1)  # silence warnings
suppressMessages({
  library(tidyverse)     # data wrangling, readr for CSV
  library(arrow)         # read_parquet
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
args.output_folder   <- args[3]  # Output dir for ALL ORF outputs
args.min_orf_length  <- if (length(args) >= 4) as.integer(args[4]) else 200L

# Create output directories
dir.create(args.output_folder, showWarnings = FALSE, recursive = TRUE)
fasta_output_dir <- file.path(args.output_folder, "fastas")
dir.create(fasta_output_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting ORF analysis...")
message("  BLAST parquet: ", args.blast_parquet)
message("  Species DB:    ", args.species_db)
message("  Output folder: ", args.output_folder)
message("  Min ORF len:   ", args.min_orf_length)

# ==============================================================================
# 2) HELPER FUNCTIONS
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

# Process a single species
process_species <- function(species, d0, args.species_db, args.output_folder,
                            fasta_output_dir, args.min_orf_length) {

  message(sprintf("Processing species: %s", species))

  # Build genome FASTA path
  genome_fasta <- file.path(args.species_db, paste0(species, ".fa"))

  # Check if genome FASTA exists
  if (!file.exists(genome_fasta)) {
    warning(sprintf("Genome FASTA not found for %s at %s. Skipping.", species, genome_fasta))
    return(NULL)
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
    return(NULL)
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

  if (is.null(genome)) return(NULL)

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
    return(NULL)
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
    warning(sprintf("No sequences with ORFs found for %s. Skipping output.", species))
    return(NULL)
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

  # Build tag vector for IDs
  tag_vec <- mcols(gr_filt)$Tag
  if (is.null(tag_vec)) tag_vec <- NA_character_
  tag_vec <- ifelse(is.na(tag_vec) | tag_vec == "", "NA", as.character(tag_vec))

  # === Export CSV ===
  csv_path <- file.path(args.output_folder, paste0(species, "_orf.csv"))
  readr::write_csv(d2_filt, csv_path)
  message(sprintf("  Wrote CSV: %s", csv_path))

  # === Export GFF3 ===
  # Drop heavy columns if present
  if ("Subject Sequence" %in% colnames(mcols(gr_filt))) {
    mcols(gr_filt)$`Subject Sequence` <- NULL
  }

  # Set feature type
  mcols(gr_filt)$type <- "region"

  # Build ID attribute
  mcols(gr_filt)$ID <- paste0(
    as.character(seqnames(gr_filt)), "_",
    start(gr_filt), "-", end(gr_filt), "|",
    tag_vec, "|",
    as.character(strand(gr_filt))
  )

  # Pass through useful BLAST columns if they exist
  keep_attr <- intersect(
    c("Query ID", "E-value", "Bit Score", "Pct Identity", "Species", "Tag"),
    colnames(mcols(gr))
  )
  if (length(keep_attr) > 0) {
    mcols(gr_filt)[keep_attr] <- mcols(gr)[has_idx, keep_attr, drop = FALSE]
  }

  gff3_path <- file.path(args.output_folder, paste0(species, "_orf.gff3"))
  rtracklayer::export(gr_filt, con = gff3_path, format = "GFF3")
  message(sprintf("  Wrote GFF3: %s", gff3_path))

  # === Export FASTA ===
  fa <- getSeq(genome, gr_filt)

  # FASTA headers: <seqid>_<start>-<end>|<tag>|<strand>
  hdr <- paste0(
    as.character(seqnames(gr_filt)), "_",
    start(gr_filt), "-", end(gr_filt), "|",
    tag_vec, "|",
    as.character(strand(gr_filt))
  )
  names(fa) <- hdr

  fasta_path <- file.path(fasta_output_dir, paste0(species, ".fa"))
  writeXStringSet(fa, fasta_path)
  message(sprintf("  Wrote FASTA: %s", fasta_path))

  return(list(
    species = species,
    total_hits = nrow(d1),
    hits_with_orfs = length(has_idx)
  ))
}

# ==============================================================================
# 3) LOAD BLAST TABLE AND GET SPECIES LIST
# ==============================================================================
message("Loading BLAST data...")
d0 <- arrow::read_parquet(args.blast_parquet)

species_list <- unique(d0$Species)
message(sprintf("Found %d species in BLAST parquet", length(species_list)))

# ==============================================================================
# 4) PROCESS ALL SPECIES
# ==============================================================================
results <- lapply(species_list, function(species) {
  tryCatch(
    process_species(species, d0, args.species_db, args.output_folder,
                    fasta_output_dir, args.min_orf_length),
    error = function(e) {
      warning(sprintf("Error processing %s: %s", species, e$message))
      return(NULL)
    }
  )
})

# Filter out NULL results
results <- Filter(Negate(is.null), results)

# ==============================================================================
# 5) SUMMARY
# ==============================================================================
if (length(results) > 0) {
  summary_df <- bind_rows(results)
  message("\n=== ORF Analysis Summary ===")
  message(sprintf("Species processed: %d", nrow(summary_df)))
  message(sprintf("Total hits analyzed: %d", sum(summary_df$total_hits)))
  message(sprintf("Total hits with ORFs: %d", sum(summary_df$hits_with_orfs)))
} else {
  message("\nNo species produced output. Check warnings above.")
}

# ==============================================================================
# 6) WRITE COMBINED PARQUET + MANIFEST
# ==============================================================================
if (length(results) > 0) {
  # Combine all per-species data into one parquet
  all_species_data <- bind_rows(lapply(species_list, function(sp) {
    csv_path <- file.path(args.output_folder, paste0(sp, "_orf.csv"))
    if (file.exists(csv_path)) {
      readr::read_csv(csv_path, show_col_types = FALSE)
    } else {
      NULL
    }
  }))

  # Write combined parquet
  arrow::write_parquet(all_species_data,
                       file.path(args.output_folder, "orf.parquet"))
  message("Wrote combined orf.parquet")

  # Write manifest (simple text file with species names that produced output)
  species_with_output <- sapply(results, function(x) x$species)
  writeLines(species_with_output,
             file.path(args.output_folder, "species_manifest.txt"))
  message("Wrote species_manifest.txt")
}

message("\nDone.")
