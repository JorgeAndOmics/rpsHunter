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
n_queries <- length(query_accessions)

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
  group_by(Species_Name, Domain_Family, Query_Accession) %>%
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
  facet_wrap(~ Query_Accession, ncol = max(1L, as.integer(n_queries))) +
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

ggsave(
  out_completeness_bars, p3,
  width     = max(12, length(domain_order_p3) * 0.55),
  height    = 4 * n_species_p3,
  dpi       = 300,
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
    )
  ) %>%
  filter(!is.na(Concordance_Quartile))

p4 <- ggplot(
  p4_data,
  aes(x = Species_Name, y = N_CDD_Domains, fill = Concordance_Quartile)
) +
  geom_violin(
    position = position_dodge(width = 0.9),
    scale    = "count",
    alpha    = 0.75
  ) +
  geom_boxplot(
    aes(group = interaction(Species_Name, Concordance_Quartile)),
    position     = position_dodge(width = 0.9),
    width        = 0.08,
    colour       = "black",
    outlier.size = 0.3
  ) +
  scale_fill_manual(values = p4_quartile_cols, name = "Concordance\nQuartile") +
  facet_wrap(~ Query_Accession, ncol = max(1L, as.integer(n_queries))) +
  labs(
    title    = "Per-Sequence Domain Complexity",
    subtitle = "Number of CDD domains per BLAST hit, stratified by concordance rate quartile",
    x = NULL,
    y = "N CDD Domains"
  ) +
  theme_rpsHunter() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold.italic"))

ggsave(out_sequence_complexity, p4, width = 14, height = 7, dpi = 300)
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
    # Extract first run of digits from chromosome name (e.g. "CM040296.1" → 40296,
    # "NC_000001.11" → 1). Using regmatches+regexpr avoids the perl=TRUE requirement
    # of gsub backreferences.
    Chr_num = suppressWarnings(
      as.integer(regmatches(Chromosome, regexpr("[0-9]+", Chromosome)))
    )
  ) %>%
  filter(!is.na(Chr_num))

p5_tile_data <- df_domains_chr %>%
  group_by(Species_Name, Chromosome, Chr_num) %>%
  summarise(
    n_domains       = n(),
    Dominant_Domain = {
      tbl <- table(Domain)
      names(tbl)[which.max(tbl)]
    },
    .groups = "drop"
  ) %>%
  mutate(log_density = log10(n_domains + 1))

# Build a globally sorted chromosome factor for a clean x-axis
chr_order <- p5_tile_data %>%
  distinct(Chromosome, Chr_num) %>%
  arrange(Chr_num, Chromosome) %>%
  pull(Chromosome)
p5_tile_data <- p5_tile_data %>%
  mutate(Chromosome = factor(Chromosome, levels = unique(chr_order)))

p5 <- ggplot(p5_tile_data, aes(x = Chromosome, y = Species_Name, fill = log_density)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(
    data = filter(p5_tile_data, n_domains > 50),
    aes(label = substr(Dominant_Domain, 1, 4)),
    size = 2.2, colour = "white"
  ) +
  scale_fill_gradient(
    low  = "#EFF3FF",
    high = "#08306B",
    name = "log10(n+1)"
  ) +
  labs(
    title    = "Chromosomal Domain Density Heatmap",
    subtitle = "log10(domain count + 1) per chromosome; dominant domain labelled where count > 50; scaffolds excluded",
    x = "Chromosome",
    y = NULL
  ) +
  theme_rpsHunter() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(face = "bold.italic")
  )

ggsave(out_chromosomal_density, p5, width = 20, height = 6, dpi = 300)
cat("  Saved:", out_chromosomal_density, "\n")

# =============================================================================
# PLOT 6 — Cross-Query Domain Comparison
# =============================================================================
cat("Plot 6: Cross-Query Domain Comparison\n")

if (n_queries >= 2) {
  acc1 <- query_accessions[1]
  acc2 <- query_accessions[2]
  lab1 <- query_labels[1]
  lab2 <- query_labels[2]

  p6_data <- df_domains %>%
    mutate(
      Domain_Family = group_domain_family(Domain),
      Query_Origin  = dplyr::case_when(
        !grepl(acc1, Query_Accession, fixed = TRUE) &
          grepl(acc2, Query_Accession, fixed = TRUE) ~ paste0(lab2, " only"),
        grepl(acc1, Query_Accession, fixed = TRUE) &
          !grepl(acc2, Query_Accession, fixed = TRUE) ~ paste0(lab1, " only"),
        grepl(acc1, Query_Accession, fixed = TRUE) &
          grepl(acc2, Query_Accession, fixed = TRUE) ~ "Both",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Query_Origin)) %>%
    group_by(Species_Name, Domain_Family, Query_Origin) %>%
    summarise(Count = n(), .groups = "drop")

  origin_levels <- c(paste0(lab1, " only"), "Both", paste0(lab2, " only"))
  p6_cols       <- c("#C0392B", "#8E44AD", "#2980B9")
  names(p6_cols) <- origin_levels

  p6_data <- p6_data %>%
    mutate(Query_Origin = factor(Query_Origin, levels = origin_levels))

  p6 <- ggplot(p6_data, aes(x = Species_Name, y = Count, fill = Query_Origin)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
    scale_fill_manual(values = p6_cols, name = "Query Origin") +
    facet_wrap(~ Domain_Family, scales = "free_y", ncol = 3) +
    labs(
      title    = "Cross-Query Domain Comparison",
      subtitle = "Domain annotations per species classified by query origin (exclusive vs shared)",
      x = NULL,
      y = "Count"
    ) +
    theme_rpsHunter() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold.italic"))

  ggsave(out_cross_query, p6, width = 16, height = 12, dpi = 300)
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

ggsave(
  out_hit_type, p7,
  width     = max(12, length(domain_order_p7) * 0.55),
  height    = 4 * n_species_p7,
  dpi       = 300,
  limitsize = FALSE
)
cat("  Saved:", out_hit_type, "\n")

cat("All 7 extended plots complete.\n")
