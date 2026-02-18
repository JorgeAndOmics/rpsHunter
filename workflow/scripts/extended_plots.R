# =============================================================================
# DEPENDENCIES
# =============================================================================
options(warn = -1)
suppressMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(yaml)
})

# =============================================================================
# INPUT / OUTPUT (injected by Snakemake script: directive)
# =============================================================================
domains_path     <- snakemake@input[["domains"]]
conc_dom_path    <- snakemake@input[["conc_domains"]]
conc_seq_path    <- snakemake@input[["conc_sequences"]]
contingency_path <- snakemake@input[["contingency"]]
config_path      <- snakemake@input[["config_path"]]

out_concordance_heatmap <- snakemake@output[["concordance_heatmap"]]
out_evidence_quality    <- snakemake@output[["evidence_quality"]]
out_completeness_bars   <- snakemake@output[["completeness_bars"]]
out_sequence_complexity <- snakemake@output[["sequence_complexity"]]
out_chromosomal_density <- snakemake@output[["chromosomal_density"]]
out_cross_query         <- snakemake@output[["cross_query"]]
out_hit_type            <- snakemake@output[["hit_type"]]

# =============================================================================
# CONFIG
# =============================================================================
config <- yaml::read_yaml(config_path)

# Species name mapping (key -> display name)
species_raw <- config$species
if (is.list(species_raw) && !is.null(names(species_raw))) {
  species_map <- unlist(species_raw)
} else {
  species_map <- setNames(unlist(species_raw), unlist(species_raw))
}

# Query accessions (accession -> label)
queries_raw <- config$queries
if (is.list(queries_raw) && !is.null(names(queries_raw))) {
  query_accessions <- names(queries_raw)
  query_labels     <- unlist(queries_raw)
} else {
  # Single-query fallback: use the single accession from config
  single_acc       <- if (!is.null(config$query)) config$query else ""
  query_accessions <- single_acc
  query_labels     <- single_acc
}
n_queries   <- length(query_accessions)
acc_to_label <- setNames(query_labels, query_accessions)   # accession → display label

# =============================================================================
# DATA LOADING (once, shared across all plots)
# =============================================================================
cat("Loading data...\n")
df_domains     <- arrow::read_parquet(domains_path)
df_conc_dom    <- arrow::read_parquet(conc_dom_path)
df_conc_seq    <- arrow::read_parquet(conc_seq_path)
df_contingency <- arrow::read_parquet(contingency_path)

cat(sprintf("  domains:               %d rows\n", nrow(df_domains)))
cat(sprintf("  concordance domains:   %d rows\n", nrow(df_conc_dom)))
cat(sprintf("  concordance sequences: %d rows\n", nrow(df_conc_seq)))
cat(sprintf("  contingency table:     %d rows\n", nrow(df_contingency)))

# =============================================================================
# SHARED HELPERS
# =============================================================================

# Custom theme
theme_rpsHunter <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 10, colour = "grey40"),
      axis.title      = element_text(face = "bold"),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.title    = element_text(face = "bold"),
      strip.text      = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# Domain family grouping via prefix matching (mirrors concordance.py logic).
# Covers all families observed in production data, including ZnF_C2H2/zf-H2C2
# zinc-finger variants that differ from the base zf-C2H2 prefix.
group_domain_family <- function(domain_vec) {
  dplyr::case_when(
    grepl("^KRAB",     domain_vec) ~ "KRAB",
    grepl("^SET",      domain_vec) ~ "SET",
    grepl("^SSXRD",    domain_vec) ~ "SSXRD",
    grepl("^zf-C2H2",  domain_vec) ~ "zf-C2H2",
    grepl("^zf-H2C2",  domain_vec) ~ "zf-H2C2",
    grepl("^ZnF_C2H2", domain_vec) ~ "ZnF_C2H2",
    grepl("^COG5048",  domain_vec) ~ "COG5048",
    grepl("^DUF4371",  domain_vec) ~ "DUF4371",
    TRUE ~ "Other"
  )
}

# =============================================================================
# PLOT 1 — Concordance Validation Heatmap
# =============================================================================
cat("Plot 1: Concordance Validation Heatmap\n")

p1_data <- df_conc_dom %>%
  filter(HMMER_Checkable == TRUE) %>%
  mutate(
    Query_Label = ifelse(
      Query_Accession %in% names(acc_to_label),
      acc_to_label[Query_Accession],
      Query_Accession
    )
  ) %>%
  group_by(Species_Name, Domain_Family, Query_Label) %>%
  summarise(
    Confirmed = sum(Concordance == "confirmed", na.rm = TRUE),
    Total     = n(),
    Rate      = Confirmed / Total,
    .groups   = "drop"
  )

p1 <- ggplot(p1_data, aes(x = Species_Name, y = Domain_Family, fill = Rate)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label = paste0(Confirmed, "/", Total)),
    size = 3, colour = "black"
  ) +
  scale_fill_viridis_c(
    name   = "Concordance\nRate",
    limits = c(0, 1),
    option = "viridis"
  ) +
  facet_wrap(~ Query_Label, ncol = max(1L, as.integer(n_queries))) +
  labs(
    title    = "Concordance Validation Heatmap",
    subtitle = "Fraction of CDD annotations confirmed by HMMER (HMMER-checkable families only)",
    x = NULL,
    y = "Domain Family"
  ) +
  theme_rpsHunter() +
  theme(axis.text.y = element_text(angle = 0))

ggsave(out_concordance_heatmap, p1, width = 14, height = 6, dpi = 300)
cat("  Saved:", out_concordance_heatmap, "\n")

# =============================================================================
# PLOT 2 — Evidence Quality Landscape
# =============================================================================
cat("Plot 2: Evidence Quality Landscape\n")

p2_data <- df_conc_dom %>%
  mutate(
    Bitscore       = suppressWarnings(as.numeric(as.character(Bitscore))),
    Blast_Bitscore = suppressWarnings(as.numeric(as.character(Blast_Bitscore)))
  ) %>%
  filter(!is.na(Bitscore) & !is.na(Blast_Bitscore) & Blast_Bitscore > 0) %>%
  filter(HMMER_Checkable == TRUE) %>%
  mutate(
    # Domain_Family already present in conc_domains (assigned by concordance.py)
    Concordance = factor(
      Concordance,
      levels = c("confirmed", "hmmer_unmatched", "hmmer_not_searched")
    )
  )

p2 <- ggplot(p2_data, aes(x = Blast_Bitscore, y = Bitscore, colour = Concordance)) +
  geom_point(alpha = 0.2, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, formula = y ~ x) +
  scale_x_log10() +
  scale_colour_manual(
    name   = "Concordance",
    values = c(
      "confirmed"          = "#1B9E77",
      "hmmer_unmatched"    = "#D95F02",
      "hmmer_not_searched" = "#AAAAAA"
    ),
    labels = c("Confirmed", "HMMER Unmatched", "Not Searched")
  ) +
  facet_wrap(~ Domain_Family, scales = "free", ncol = 4) +
  labs(
    title    = "Evidence Quality Landscape",
    subtitle = "CDD bitscore vs. BLAST bitscore, coloured by HMMER concordance",
    x = "BLAST Bitscore (log10 scale)",
    y = "CDD Bitscore"
  ) +
  theme_rpsHunter()

ggsave(out_evidence_quality, p2, width = 16, height = 10, dpi = 300)
cat("  Saved:", out_evidence_quality, "\n")

# =============================================================================
# PLOT 3 — Domain Completeness Proportional Bar Chart
# =============================================================================
cat("Plot 3: Domain Completeness Proportional Bars\n")

# contingency_table.parquet already carries mapped Incomplete labels
p3_fill_cols <- c(
  "Complete"    = "#2166AC",
  "N-Truncated" = "#4DAC26",
  "C-Truncated" = "#F4A582",
  "Bitruncated" = "#D6604D"
)

# Order domains by total count descending
domain_order_p3 <- df_contingency %>%
  group_by(Domain) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Domain)

p3_data <- df_contingency %>%
  filter(Count > 0) %>%
  group_by(Species_Name, Domain) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  mutate(
    Domain     = factor(Domain, levels = domain_order_p3),
    Incomplete = factor(
      Incomplete,
      levels = c("Complete", "N-Truncated", "C-Truncated", "Bitruncated")
    ),
    Bar_Label = ifelse(Proportion > 0.05, as.character(as.integer(Count)), "")
  )

# PNG devices hard-cap at 50,000 px per dimension regardless of limitsize = FALSE.
# For plots that scale height with species count we compute DPI dynamically so
# height_in * dpi stays within 49,000 px, with a floor of 72 DPI.
MAX_PNG_PX  <- 49000L

n_species_p3 <- length(unique(p3_data$Species_Name))

p3 <- ggplot(p3_data, aes(x = Domain, y = Count, fill = Incomplete)) +
  geom_col(position = "fill", width = 0.8) +
  geom_text(
    aes(label = Bar_Label),
    position = position_fill(vjust = 0.5),
    size = 2.5, colour = "white"
  ) +
  scale_fill_manual(values = p3_fill_cols, name = "Completeness", drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~ Species_Name, ncol = 1, scales = "free_x") +
  labs(
    title    = "Domain Completeness Proportional Bar Chart",
    subtitle = "Proportion of Complete / Truncated / Bitruncated annotations per domain, per species",
    x = "Domain",
    y = "Proportion"
  ) +
  theme_rpsHunter() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    strip.text  = element_text(face = "bold.italic", size = 9)
  )

p3_height_in <- 4 * n_species_p3
p3_dpi       <- max(72L, min(300L, as.integer(MAX_PNG_PX / p3_height_in)))
ggsave(
  out_completeness_bars, p3,
  width     = max(12, length(domain_order_p3) * 0.55),
  height    = p3_height_in,
  dpi       = p3_dpi,
  limitsize = FALSE
)
cat("  Saved:", out_completeness_bars, "\n")

# =============================================================================
# PLOT 4 — Per-Sequence Domain Complexity
# =============================================================================
cat("Plot 4: Per-Sequence Domain Complexity\n")

p4_quartile_cols <- c(
  "Q1 (0-25%)"   = "#D73027",
  "Q2 (25-50%)"  = "#FC8D59",
  "Q3 (50-75%)"  = "#91CF60",
  "Q4 (75-100%)" = "#1A9850"
)

p4_data <- df_conc_seq %>%
  filter(N_CDD_Domains > 0 & !is.na(Concordance_Rate)) %>%
  mutate(
    Concordance_Quartile = cut(
      Concordance_Rate,
      breaks         = c(-Inf, 0.25, 0.5, 0.75, Inf),
      labels         = c("Q1 (0-25%)", "Q2 (25-50%)", "Q3 (50-75%)", "Q4 (75-100%)"),
      right          = TRUE,
      include.lowest = TRUE
    ),
    Query_Label = ifelse(
      Query_Accession %in% names(acc_to_label),
      acc_to_label[Query_Accession],
      Query_Accession
    )
  ) %>%
  filter(!is.na(Concordance_Quartile))

# Design: x = quartile (4 clearly-labelled categories), facet_grid rows=species / cols=query.
# This gives one violin per quartile per panel — every x-axis tick maps to exactly one violin.
p4 <- ggplot(
  p4_data,
  aes(x = Concordance_Quartile, y = N_CDD_Domains, fill = Concordance_Quartile)
) +
  geom_violin(scale = "count", alpha = 0.75, trim = TRUE) +
  geom_boxplot(
    width        = 0.10,
    colour       = "black",
    fill         = "white",
    outlier.size = 0.3,
    outlier.alpha = 0.4
  ) +
  scale_fill_manual(values = p4_quartile_cols, name = "Concordance\nQuartile", guide = "none") +
  facet_grid(rows = vars(Species_Name), cols = vars(Query_Label)) +
  labs(
    title    = "Per-Sequence Domain Complexity by Concordance Quartile",
    subtitle = paste0(
      "Each panel = one species \u00d7 one query. ",
      "x-axis = concordance rate quartile of the BLAST hit sequence. ",
      "y-axis = number of CDD domains annotated on that sequence."
    ),
    x = "Concordance Rate Quartile",
    y = "N CDD Domains per Sequence"
  ) +
  theme_rpsHunter() +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1),
    strip.text.y = element_text(face = "bold.italic", angle = 0, hjust = 0),
    strip.text.x = element_text(face = "bold")
  )

p4_height_in <- 3 * length(unique(p4_data$Species_Name))
p4_dpi       <- max(72L, min(300L, as.integer(MAX_PNG_PX / p4_height_in)))
ggsave(out_sequence_complexity, p4, width = 10, height = p4_height_in, dpi = p4_dpi, limitsize = FALSE)
cat("  Saved:", out_sequence_complexity, "\n")

# =============================================================================
# PLOT 5 — Chromosomal Domain Density Heatmap
# =============================================================================
cat("Plot 5: Chromosomal Domain Density Heatmap\n")

# Ensure Species_Name is present
if (!"Species_Name" %in% names(df_domains)) {
  df_domains <- df_domains %>%
    mutate(
      Species_Name = ifelse(
        Species %in% names(species_map),
        species_map[Species],
        Species
      )
    )
}

df_domains_chr <- df_domains %>%
  filter(!is.na(Chromosome) & nzchar(Chromosome)) %>%
  filter(grepl("[0-9]", Chromosome)) %>%
  mutate(
    # Extract the first long run of digits for within-species ordering.
    # For CM accessions (e.g. "CM040296.1") this gives the assembly-sequential
    # integer; for NC accessions (e.g. "NC_000001.11") it gives 1.
    Chr_num = suppressWarnings(
      as.integer(regmatches(Chromosome, regexpr("[0-9]+", Chromosome)))
    ),
    # Group domains into families using the same logic as group_domain_family()
    Domain_Family = dplyr::case_when(
      grepl("^KRAB",     Domain) ~ "KRAB",
      grepl("^SET",      Domain) ~ "SET",
      grepl("^SSXRD",    Domain) ~ "SSXRD",
      grepl("^zf-C2H2",  Domain) ~ "zf-C2H2",
      grepl("^zf-H2C2",  Domain) ~ "zf-H2C2",
      grepl("^ZnF_C2H2", Domain) ~ "ZnF_C2H2",
      grepl("^COG5048",  Domain) ~ "COG5048",
      grepl("^DUF4371",  Domain) ~ "DUF4371",
      TRUE ~ "Other"
    )
  ) %>%
  filter(!is.na(Chr_num))

# Count per species × chromosome × domain family
p5_stacked_data <- df_domains_chr %>%
  group_by(Species_Name, Chromosome, Chr_num, Domain_Family) %>%
  summarise(Count = n(), .groups = "drop")

# Per-species chromosome ordering by Chr_num (assembly position)
p5_stacked_data <- p5_stacked_data %>%
  group_by(Species_Name) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = unique(Chromosome[order(Chr_num)]))
  ) %>%
  ungroup()

# Consistent domain family factor order (by total abundance, most frequent first)
family_order_p5 <- p5_stacked_data %>%
  group_by(Domain_Family) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Domain_Family)
p5_stacked_data <- p5_stacked_data %>%
  mutate(Domain_Family = factor(Domain_Family, levels = family_order_p5))

domain_family_cols <- c(
  "zf-C2H2"   = "#984EA3",
  "COG5048"   = "#377EB8",
  "SET"       = "#4DAF4A",
  "zf-H2C2"   = "#FF7F00",
  "ZnF_C2H2"  = "#A65628",
  "KRAB"      = "#E41A1C",
  "SSXRD"     = "#F781BF",
  "DUF4371"   = "#999999",
  "Other"     = "#DDDDDD"
)

n_species_p5 <- length(unique(p5_stacked_data$Species_Name))

p5 <- ggplot(p5_stacked_data, aes(x = Chromosome, y = Count, fill = Domain_Family)) +
  geom_col(position = "stack", width = 0.9) +
  scale_fill_manual(
    values = domain_family_cols,
    name   = "Domain Family",
    drop   = FALSE
  ) +
  facet_wrap(~ Species_Name, ncol = 1, scales = "free_x") +
  labs(
    title    = "Chromosomal Domain Distribution",
    subtitle = paste0(
      "Stacked count of all domain families per chromosome per species. ",
      "Chromosomes ordered by assembly position (accession order). Scaffolds excluded."
    ),
    x = "Chromosome (assembly accession)",
    y = "Domain Count"
  ) +
  theme_rpsHunter() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    strip.text  = element_text(face = "bold.italic", size = 10)
  )

p5_height_in <- 4 * n_species_p5
p5_dpi       <- max(72L, min(300L, as.integer(MAX_PNG_PX / p5_height_in)))
ggsave(
  out_chromosomal_density, p5,
  width     = 20,
  height    = p5_height_in,
  dpi       = p5_dpi,
  limitsize = FALSE
)
cat("  Saved:", out_chromosomal_density, "\n")

# =============================================================================
# PLOT 6 — Cross-Query Domain Comparison
# =============================================================================
cat("Plot 6: Cross-Query Domain Comparison\n")

if (n_queries >= 2) {
  # For each domain row, detect which of the N configured query accessions appear
  # in the (possibly comma-separated) Query_Accession field, map to labels, then
  # classify as exclusive (one query) or shared (>1 queries).
  p6_data <- df_domains %>%
    mutate(Domain_Family = group_domain_family(Domain)) %>%
    rowwise() %>%
    mutate(
      present_labels = list(
        query_labels[vapply(query_accessions,
                            function(acc) grepl(acc, Query_Accession, fixed = TRUE),
                            logical(1))]
      ),
      N_Queries_Hit = length(present_labels),
      Query_Origin  = dplyr::case_when(
        N_Queries_Hit == 0 ~ NA_character_,
        N_Queries_Hit == 1 ~ present_labels[[1]],
        TRUE               ~ paste0("Shared (", N_Queries_Hit, " queries)")
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(Query_Origin)) %>%
    group_by(Species_Name, Domain_Family, Query_Origin) %>%
    summarise(Count = n(), .groups = "drop")

  # Colour palette: one colour per query label + grey scale for shared categories
  query_palette <- c("#C0392B", "#2980B9", "#27AE60", "#8E44AD",
                     "#E67E22", "#16A085", "#2C3E50", "#8E44AD")
  exclusive_cols <- setNames(query_palette[seq_len(n_queries)], query_labels)

  shared_levels <- sort(unique(grep("^Shared", p6_data$Query_Origin, value = TRUE)))
  n_shared      <- length(shared_levels)
  shared_cols   <- setNames(
    grey.colors(max(1L, n_shared), start = 0.55, end = 0.80),
    shared_levels
  )

  p6_cols       <- c(exclusive_cols, shared_cols)
  origin_levels <- c(query_labels, shared_levels)

  p6_data <- p6_data %>%
    mutate(Query_Origin = factor(Query_Origin, levels = origin_levels))

  p6 <- ggplot(p6_data, aes(x = Species_Name, y = Count, fill = Query_Origin)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
    scale_fill_manual(values = p6_cols, name = "Query Origin") +
    facet_wrap(~ Domain_Family, scales = "free_y", ncol = 3) +
    labs(
      title    = "Cross-Query Domain Comparison",
      subtitle = paste0(
        "Domain annotations per species classified by query origin. ",
        "Solid colours = exclusive to one query; grey = detected by multiple queries."
      ),
      x = NULL,
      y = "Count"
    ) +
    theme_rpsHunter() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold.italic"))

  ggsave(out_cross_query, p6, width = 16, height = 12, dpi = 300, limitsize = FALSE)
  cat("  Saved:", out_cross_query, "\n")

} else {
  # Single-query fallback
  p6_blank <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label  = "Cross-query comparison requires\nat least 2 configured queries.",
      size   = 7, hjust = 0.5, vjust = 0.5, colour = "grey50"
    ) +
    theme_void() +
    labs(title = "Cross-Query Domain Comparison (N/A — single query)")
  ggsave(out_cross_query, p6_blank, width = 10, height = 6, dpi = 300)
  cat("  Saved (blank — single query):", out_cross_query, "\n")
}

# =============================================================================
# PLOT 7 — Hit-Type Distribution by Domain
# =============================================================================
cat("Plot 7: Hit-Type Distribution by Domain\n")

p7_data <- df_domains %>%
  filter(!is.na(Hit_type) & nzchar(as.character(Hit_type))) %>%
  mutate(
    Hit_type = factor(
      Hit_type,
      levels = c("Specific", "Non-specific", "Superfamily")
    )
  ) %>%
  group_by(Species_Name, Domain, Hit_type) %>%
  summarise(Count = n(), .groups = "drop")

# Order Domain by total Specific-hit count descending; non-Specific domains append at end
domain_specific_order <- p7_data %>%
  filter(Hit_type == "Specific") %>%
  group_by(Domain) %>%
  summarise(Specific_Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Specific_Total)) %>%
  pull(Domain)

all_domains_p7   <- unique(p7_data$Domain)
remaining_p7     <- setdiff(all_domains_p7, domain_specific_order)
domain_order_p7  <- c(domain_specific_order, remaining_p7)

p7_data <- p7_data %>%
  mutate(Domain = factor(Domain, levels = domain_order_p7))

p7_fill_cols <- c(
  "Specific"     = "#1F3B6E",
  "Non-specific" = "#5B8DB8",
  "Superfamily"  = "#AECDE0"
)

n_species_p7 <- length(unique(p7_data$Species_Name))

p7 <- ggplot(p7_data, aes(x = Domain, y = Count, fill = Hit_type)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = p7_fill_cols, name = "Hit Type", drop = FALSE) +
  facet_wrap(~ Species_Name, ncol = 1, scales = "free_y") +
  labs(
    title    = "Hit-Type Distribution by Domain",
    subtitle = "Stacked count of Specific / Non-specific / Superfamily CDD hits per domain per species",
    x = "Domain",
    y = "Count"
  ) +
  theme_rpsHunter() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    strip.text  = element_text(face = "bold.italic", size = 9)
  )

p7_height_in <- 4 * n_species_p7
p7_dpi       <- max(72L, min(300L, as.integer(MAX_PNG_PX / p7_height_in)))
ggsave(
  out_hit_type, p7,
  width     = max(12, length(domain_order_p7) * 0.55),
  height    = p7_height_in,
  dpi       = p7_dpi,
  limitsize = FALSE
)
cat("  Saved:", out_hit_type, "\n")

cat("All 7 extended plots complete.\n")
