#!/usr/bin/env Rscript
# =============================================================================
# merge_queries.R — Cross-query domain deduplication via GenomicRanges
#
# Usage:
#   Rscript merge_queries.R <output.parquet> <config.yaml> <input1.parquet> [input2.parquet ...]
#
# Loads all per-query domain parquets, concatenates, groups by
# Species × Chromosome × Domain, merges overlapping [Start, End] intervals
# via GenomicRanges reduce(with.revmap=TRUE), and applies collapse rules.
# =============================================================================

options(warn = -1)
suppressMessages({
  library(arrow)
  library(tidyverse)
  library(GenomicRanges)
  library(IRanges)
  library(yaml)
})

# =============================================================================
# Command-line arguments
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: merge_queries.R <output.parquet> <config.yaml> <input1.parquet> [input2.parquet ...]")
}

output_parquet <- args[1]
config_path    <- args[2]
input_files    <- args[3:length(args)]

# =============================================================================
# Load and concatenate per-query domain tables
# =============================================================================
cat("Loading per-query domain tables...\n")

dfs <- list()
for (f in input_files) {
  if (file.exists(f)) {
    df <- arrow::read_parquet(f)
    if (nrow(df) > 0) {
      dfs <- c(dfs, list(df))
      cat(sprintf("  %s: %d rows\n", basename(f), nrow(df)))
    }
  }
}

if (length(dfs) == 0) {
  cat("No input data — writing empty parquet.\n")
  dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(tibble(), output_parquet)
  # Also write CSV companion
  write_csv(tibble(), sub("\\.parquet$", ".csv", output_parquet))
  quit(save = "no", status = 0)
}

combined <- bind_rows(dfs)
cat(sprintf("Combined: %d rows from %d files\n", nrow(combined), length(dfs)))

# If only one query, skip deduplication — just pass through
unique_queries <- unique(combined$Query_Accession)
if (length(unique_queries) <= 1) {
  cat("Single query detected — skipping deduplication, writing as-is.\n")
  dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(combined, output_parquet)
  write_csv(combined, sub("\\.parquet$", ".csv", output_parquet))
  cat(sprintf("Output: %d rows\n", nrow(combined)))
  quit(save = "no", status = 0)
}

# =============================================================================
# Prepare data for GenomicRanges
# =============================================================================
cat("Preparing genomic ranges for deduplication...\n")

# Clean numeric columns and filter out rows without valid coordinates
combined <- combined %>%
  mutate(
    Start_num = as.integer(str_extract(as.character(Start), "^[0-9]+")),
    End_num   = as.integer(str_extract(as.character(End), "^[0-9]+")),
    Bitscore_num = as.numeric(Bitscore),
    Evalue_num   = as.numeric(Evalue),
    From_num     = as.integer(From),
    To_num       = as.integer(To),
    Seq_length_num = as.integer(Seq_length),
    .row_id = row_number()
  )

# Rows with valid coordinates for range-based dedup
has_coords <- !is.na(combined$Start_num) & !is.na(combined$End_num) &
              nzchar(combined$Chromosome) & !is.na(combined$Chromosome) &
              nzchar(combined$Domain) & !is.na(combined$Domain)

dedup_df   <- combined %>% filter(has_coords)
nocoord_df <- combined %>% filter(!has_coords)

if (nrow(dedup_df) == 0) {
  cat("No rows with valid coordinates — writing combined as-is.\n")
  out <- combined %>% select(-.row_id, -Start_num, -End_num, -Bitscore_num,
                              -Evalue_num, -From_num, -To_num, -Seq_length_num)
  dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(out, output_parquet)
  write_csv(out, sub("\\.parquet$", ".csv", output_parquet))
  cat(sprintf("Output: %d rows\n", nrow(out)))
  quit(save = "no", status = 0)
}

# Fix start/end directionality
dedup_df <- dedup_df %>%
  mutate(
    Range_Start = pmin(Start_num, End_num),
    Range_End   = pmax(Start_num, End_num)
  )

# =============================================================================
# Hit_type and Incomplete priority maps
# =============================================================================
hit_type_priority <- c("Specific" = 1, "Non-specific" = 2, "Superfamily" = 3)
incomplete_priority <- c("-" = 1, "C" = 2, "N" = 2, "NC" = 3)

best_hit_type <- function(vals) {
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) return("")
  ranks <- hit_type_priority[vals]
  ranks[is.na(ranks)] <- 99
  vals[which.min(ranks)]
}

best_incomplete <- function(vals) {
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) return("")
  ranks <- incomplete_priority[vals]
  ranks[is.na(ranks)] <- 99
  vals[which.min(ranks)]
}

paste_unique <- function(vals) {
  vals <- vals[!is.na(vals) & nzchar(vals)]
  paste(unique(vals), collapse = ",")
}

# =============================================================================
# Group by Species × Chromosome × Domain, reduce overlapping ranges
# =============================================================================
cat("Building GenomicRanges and reducing overlaps...\n")

# Create a grouping key
dedup_df <- dedup_df %>%
  mutate(group_key = paste(Species, Chromosome, Domain, sep = "|||"))

group_keys <- unique(dedup_df$group_key)
cat(sprintf("  %d unique Species x Chromosome x Domain groups\n", length(group_keys)))

result_rows <- list()

for (gk in group_keys) {
  grp <- dedup_df %>% filter(group_key == gk)

  # Build GRanges for this group
  gr <- GRanges(
    seqnames = grp$Chromosome,
    ranges   = IRanges(start = grp$Range_Start, end = grp$Range_End)
  )
  mcols(gr)$.row_id <- grp$.row_id

  # Reduce with revmap to track which original rows contributed
  gr_reduced <- reduce(gr, with.revmap = TRUE)

  for (j in seq_along(gr_reduced)) {
    contributing_idx <- mcols(gr_reduced)$revmap[[j]]
    contributing_rows <- grp[contributing_idx, ]

    # Best bitscore row for Definition
    best_bs_idx <- which.max(contributing_rows$Bitscore_num)

    merged_row <- tibble(
      Species              = contributing_rows$Species[1],
      Species_Name         = contributing_rows$Species_Name[1],
      Query_Accession      = paste_unique(contributing_rows$Query_Accession),
      Session_ordinal      = paste_unique(contributing_rows$Session_ordinal),
      Program              = contributing_rows$Program[1],
      Version              = contributing_rows$Version[1],
      Database             = contributing_rows$Database[1],
      Score_matrix         = contributing_rows$Score_matrix[1],
      Evalue_threshold     = contributing_rows$Evalue_threshold[1],
      Query_ID             = paste_unique(contributing_rows$Query_ID),
      Seq_type             = contributing_rows$Seq_type[1],
      Seq_length           = as.character(max(contributing_rows$Seq_length_num, na.rm = TRUE)),
      Definition           = contributing_rows$Definition[best_bs_idx],
      Chromosome           = contributing_rows$Chromosome[1],
      Start                = as.character(start(gr_reduced[j])),
      End                  = as.character(end(gr_reduced[j])),
      Hit_type             = best_hit_type(contributing_rows$Hit_type),
      PSSM_ID              = contributing_rows$PSSM_ID[1],
      From                 = as.character(min(contributing_rows$From_num, na.rm = TRUE)),
      To                   = as.character(max(contributing_rows$To_num, na.rm = TRUE)),
      Evalue               = as.character(min(contributing_rows$Evalue_num, na.rm = TRUE)),
      Bitscore             = as.character(max(contributing_rows$Bitscore_num, na.rm = TRUE)),
      Accession            = contributing_rows$Accession[1],
      Domain               = contributing_rows$Domain[1],
      Incomplete           = best_incomplete(contributing_rows$Incomplete),
      Superfamily_PSSM_ID  = contributing_rows$Superfamily_PSSM_ID[1],
      Tag                  = paste_unique(contributing_rows$Tag)
    )

    result_rows <- c(result_rows, list(merged_row))
  }
}

# =============================================================================
# Combine results
# =============================================================================
deduped <- bind_rows(result_rows)

# Add back rows without coordinates (if any)
if (nrow(nocoord_df) > 0) {
  nocoord_out <- nocoord_df %>%
    select(-.row_id, -Start_num, -End_num, -Bitscore_num,
           -Evalue_num, -From_num, -To_num, -Seq_length_num)

  # Ensure Tag column exists in nocoord output
  if (!"Tag" %in% colnames(nocoord_out)) {
    nocoord_out$Tag <- ""
  }

  deduped <- bind_rows(deduped, nocoord_out)
}

cat(sprintf("Deduplication: %d -> %d rows (%.1f%% reduction)\n",
            nrow(combined), nrow(deduped),
            (1 - nrow(deduped) / nrow(combined)) * 100))

# =============================================================================
# Write output
# =============================================================================
dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)
arrow::write_parquet(deduped, output_parquet)
write_csv(deduped, sub("\\.parquet$", ".csv", output_parquet))
cat(sprintf("Output written: %s (%d rows)\n", output_parquet, nrow(deduped)))
