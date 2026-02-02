################################################################################
# ORF Analysis - Enriches blast.parquet with ORF detection results
################################################################################

# ==============================================================================
# SECTION 1: DEPENDENCIES
# ==============================================================================
options(warn = -1)
suppressMessages({
  library(argparse)
  library(yaml)
  library(tidyverse)
  library(arrow)
  library(GenomicRanges)
  library(Biostrings)
  library(ORFik)
  library(IRanges)
})

# ==============================================================================
# SECTION 2: ARGUMENT PARSING
# ==============================================================================
parser <- ArgumentParser(description = "Enrich blast.parquet with ORF detection results")

parser$add_argument("blast_parquet",
  help = "Path to blast.parquet file (will be modified in place)")

parser$add_argument("species_db",
  help = "Directory containing genome FASTA files ({Species}.fa)")

parser$add_argument("output_folder",
  help = "Output directory for flag file and species manifest")

args <- parser$parse_args()

# ==============================================================================
# SECTION 3: CONFIGURATION
# ==============================================================================
# Read ORF parameters from config.yaml
CONFIG_PATH <- file.path(dirname(args$blast_parquet), "..", "..", "data", "config", "config.yaml")
config <- yaml::read_yaml(CONFIG_PATH)

# ORF parameters
MIN_ORF_LENGTH <- config$orf$min_orf_length %||% 200L
START_CODON    <- {
  codon <- config$orf$start_codon %||% ""
  if (nchar(codon) == 3 && grepl("^[ATGC]{3}$", codon)) codon else NULL
}
LONGEST_ORF    <- config$orf$longest_orf %||% TRUE

# Derived paths
BLAST_DIR      <- dirname(args$blast_parquet)
BLAST_CSV      <- file.path(BLAST_DIR, "blast.csv")
FLAG_FILE      <- file.path(args$output_folder, "orf_analysis.flag")
MANIFEST_FILE  <- file.path(args$output_folder, "species_manifest.txt")

dir.create(args$output_folder, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# SECTION 4: HELPER FUNCTIONS
# ==============================================================================

#' Normalize sequence IDs for matching between BLAST and FASTA
normalize_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub(" .*$", "", x)              # Drop FASTA description
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- sub("\\.\\d+$", "", x)          # Drop version suffix
  x
}

#' Summarize ORFs as string: "seqid:start-end(strand):AAlen"
summarize_orfs <- function(orf_iranges, parent_granges) {
  if (length(orf_iranges) == 0) return(NA_character_)

  local_start <- start(orf_iranges)
  local_end   <- end(orf_iranges)
  str         <- as.character(strand(parent_granges))

  # Convert local coords to genomic coords
  if (str == "+") {
    genome_start <- start(parent_granges) + local_start - 1L
    genome_end   <- start(parent_granges) + local_end - 1L
  } else {
    genome_start <- end(parent_granges) - local_end + 1L
    genome_end   <- end(parent_granges) - local_start + 1L
  }

  aa_len <- floor((local_end - local_start + 1L) / 3L)

  paste0(seqnames(parent_granges), ":", genome_start, "-", genome_end,
         "(", str, "):", aa_len, "aa", collapse = ";")
}

# ==============================================================================
# SECTION 5: SPECIES PROCESSING
# ==============================================================================

#' Process one species: find ORFs and return Tag + ORF metadata
#' @return tibble(Tag, ORF_COUNT, ORF_SUMMARY) or NULL if no ORFs found
process_species <- function(species, blast_data, species_db) {

  # --- 5.1 Filter BLAST hits for this species ---
  hits <- blast_data %>%
    filter(Species == species, !is.na(`Subject ID`), !is.na(`S. Start`), !is.na(`S. End`)) %>%
    mutate(
      strand = if_else(`S. Start` < `S. End`, "+", "-"),
      start  = pmin(as.integer(`S. Start`), as.integer(`S. End`)),
      end    = pmax(as.integer(`S. Start`), as.integer(`S. End`))
    )

  if (nrow(hits) == 0) return(NULL)

  # --- 5.2 Load genome FASTA ---
  genome_path <- file.path(species_db, paste0(species, ".fa"))
  if (!file.exists(genome_path)) {
    warning("Genome not found: ", genome_path)
    return(NULL)
  }
  genome <- readDNAStringSet(genome_path)

  # --- 5.3 Build GRanges and harmonize IDs ---
  gr <- makeGRangesFromDataFrame(hits,
    seqnames.field = "Subject ID", start.field = "start",
    end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)

  # Map BLAST IDs to FASTA IDs via normalization
  blast_ids_norm  <- normalize_ids(seqlevels(gr))
  genome_ids_norm <- normalize_ids(names(genome))
  mapped <- names(genome)[match(blast_ids_norm, genome_ids_norm)]

  valid <- !is.na(mapped)
  if (any(valid)) {
    gr <- renameSeqlevels(gr, setNames(mapped[valid], seqlevels(gr)[valid]))
  }
  gr <- keepSeqlevels(gr, intersect(seqlevels(gr), names(genome)), pruning.mode = "coarse")

  if (length(gr) == 0) return(NULL)

  # --- 5.4 Extract sequences and find ORFs ---
  seqs <- getSeq(genome, gr)

  orf_args <- list(seqs = seqs, minimumLength = MIN_ORF_LENGTH, longestORF = LONGEST_ORF)
  if (!is.null(START_CODON)) orf_args$startCodon <- START_CODON
  orfs <- do.call(ORFik::findORFs, orf_args)

  # --- 5.5 Filter to sequences with ORFs ---
  orf_counts <- elementNROWS(orfs)
  has_orf <- which(orf_counts > 0)

  if (length(has_orf) == 0) return(NULL)

  message(sprintf("  %s: %d/%d sequences have ORFs", species, length(has_orf), length(gr)))

  # --- 5.6 Build result tibble ---
  gr_filt   <- gr[has_orf]
  orfs_filt <- orfs[has_orf]

  tibble(
    Tag = mcols(gr_filt)$Tag,
    ORF_COUNT = as.integer(orf_counts[has_orf]),
    ORF_SUMMARY = vapply(seq_along(gr_filt),
      function(i) summarize_orfs(orfs_filt[[i]], gr_filt[i]), character(1))
  )
}

# ==============================================================================
# SECTION 6: MAIN EXECUTION
# ==============================================================================

message("ORF Analysis")
message("  Config: min_length=", MIN_ORF_LENGTH,
        ", start_codon=", START_CODON %||% "any",
        ", longest_orf=", LONGEST_ORF)

# --- 6.1 Load BLAST data ---
blast_data <- read_parquet(args$blast_parquet)

# Remove previous ORF columns if re-running
if ("ORF_Filtered" %in% names(blast_data)) {
  message("  Removing previous ORF columns...")
  blast_data <- select(blast_data, -any_of(c("ORF_Filtered", "ORF_COUNT", "ORF_SUMMARY")))
}

species_list <- unique(blast_data$Species)
message("  Species: ", length(species_list))

# --- 6.2 Process all species ---
orf_results <- lapply(species_list, function(sp) {
  tryCatch(process_species(sp, blast_data, args$species_db), error = function(e) NULL)
})
orf_info <- bind_rows(Filter(Negate(is.null), orf_results))

# --- 6.3 Enrich BLAST data with ORF columns ---
if (nrow(orf_info) == 0) {
  message("  No ORFs found")
  blast_enriched <- mutate(blast_data,
    ORF_Filtered = FALSE, ORF_COUNT = NA_integer_, ORF_SUMMARY = NA_character_)
} else {
  blast_enriched <- blast_data %>%
    left_join(orf_info, by = "Tag") %>%
    mutate(ORF_Filtered = !is.na(ORF_COUNT), ORF_COUNT = as.integer(ORF_COUNT))
}

n_with_orf <- sum(blast_enriched$ORF_Filtered)
message("  Result: ", n_with_orf, "/", nrow(blast_enriched), " sequences with ORFs")

# ==============================================================================
# SECTION 7: OUTPUT
# ==============================================================================

# --- 7.1 Write enriched BLAST tables ---
write_parquet(blast_enriched, args$blast_parquet)
write_csv(blast_enriched, BLAST_CSV)

# --- 7.2 Write flag file ---
writeLines(c(
  paste0("timestamp: ", Sys.time()),
  paste0("min_orf_length: ", MIN_ORF_LENGTH),
  paste0("start_codon: ", START_CODON %||% ""),
  paste0("longest_orf: ", LONGEST_ORF),
  paste0("sequences_with_orfs: ", n_with_orf),
  paste0("status: ", if (n_with_orf > 0) "complete" else "no_orfs_found")
), FLAG_FILE)

# --- 7.3 Write species manifest ---
species_with_orfs <- unique(blast_enriched$Species[blast_enriched$ORF_Filtered])
writeLines(species_with_orfs, MANIFEST_FILE)

message("Done")
