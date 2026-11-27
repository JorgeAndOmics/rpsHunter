# The pages of the figure PDFs. Each function takes prepared data and the
# layout and returns a page (a plot with its size). Species are always on
# rows, in one order, and every configured species keeps its row.

# Agreement labels a CDD call can have (HMMER-only applies to HMMER's own table).
CALL_AGREEMENT <- c("both", "cdd_only", "not_compared")

#' The parts every species-by-family page shares: one panel per family, italic
#' species rows that are never dropped, the theme.
species_by_family <- function(plot) {
  plot +
    family_panels(scales = "free_x") +
    scale_y_discrete(drop = FALSE) +
    labs(y = NULL) +
    theme_rpshunter() +
    italic_rows()
}

#' One panel per family, titled in words, empty families kept.
family_panels <- function(scales = "fixed") {
  facet_wrap(
    vars(Family), nrow = 1, scales = scales, drop = FALSE, labeller = word_labeller()
  )
}

#' Size of a species-by-family page.
species_page <- function(plot, layout, panels = length(layout$families)) {
  page(plot, page_size(rows = length(layout$species), panels = panels))
}

#' Calls with their species and family as ordered factors.
laid_out <- function(data, layout) {
  mutate(
    data,
    Species = species_rows(Species, layout),
    Family = family_columns(Family, layout)
  )
}

# ---- overview ----------------------------------------------------------------

page_family_presence <- function(run, layout) {
  presence <- family_presence(run, layout) |>
    mutate(
      Species = species_rows(Species, layout),
      Family = family_columns(Family, layout, extra = ALL_FAMILIES)
    )
  plot <- ggplot(presence, aes(x = Family, y = Species, fill = N)) +
    geom_tile(colour = PAPER, linewidth = 1) +
    tile_numbers(size = 3.2) +
    scale_fill_amount(guide = "none") +
    scale_x_discrete(position = "top", drop = FALSE, labels = words) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = "Loci carrying each domain family",
      subtitle = paste(
        "A locus is a group of query-like hits within the chaining gap.",
        "The last column counts loci that carry every configured family."
      ),
      x = NULL, y = NULL
    ) +
    theme_rpshunter() +
    italic_rows() +
    theme(panel.grid = element_blank())
  page(plot, page_size(length(layout$species), length(layout$families) + 1, 0.3, 1.1))
}

page_funnel <- function(run, layout) {
  counts <- funnel_counts(run, layout) |>
    mutate(Species = species_rows(Species, layout))
  plot <- ggplot(counts, aes(x = N, y = Species)) +
    geom_col(fill = RAMP_DARK, width = 0.7) +
    geom_text(
      aes(label = count_labels(N)),
      hjust = -0.15, family = FONT, size = 2.9, colour = SOFT_INK
    ) +
    facet_wrap(vars(Step), nrow = 1, scales = "free_x") +
    scale_x_continuous(labels = count_labels, expand = expansion(mult = c(0, 0.35))) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = "From hits to loci",
      subtitle = paste(
        "What remains at each step. Hits are selected by quality alone;",
        "a locus chains selected hits; a family needs at least one domain call."
      ),
      x = NULL, y = NULL
    ) +
    theme_rpshunter() +
    italic_rows() +
    theme(panel.grid.major.y = element_blank())
  page(plot, page_size(length(layout$species), length(FUNNEL_STEPS), panel_width = 1.9))
}

# ---- domains -------------------------------------------------------------------

#' A species-by-family page of stacked bars, as counts or as shares.
page_stacked <- function(calls, layout, column, fill_scale, shares, title, subtitle) {
  data <- laid_out(calls, layout)
  plot <- ggplot(data, aes(y = Species, fill = .data[[column]])) +
    geom_bar(
      position = if (shares) "fill" else "stack", width = 0.7, show.legend = TRUE
    ) +
    fill_scale +
    scale_x_continuous(
      labels = if (shares) scales::label_percent() else count_labels,
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(title = title, subtitle = subtitle, x = NULL)
  species_page(species_by_family(plot), layout)
}

page_evidence <- function(run, layout) {
  calls <- mutate(run$calls, Evidence = factor(Evidence, EVIDENCE_LEVELS))
  page_stacked(
    calls, layout, "Evidence", scale_fill_evidence(name = "Evidence tier"),
    shares = FALSE,
    title = "Domain calls by strength of evidence",
    subtitle = paste(
      "CDD calls per family. E-values are on the scale of the full CDD,",
      "as NCBI's CD-Search reports them."
    )
  )
}

page_hit_type <- function(run, layout) {
  calls <- mutate(run$calls, Hit_Type = factor(Hit_Type, HIT_TYPE_LEVELS))
  page_stacked(
    calls, layout, "Hit_Type", scale_fill_hit_type(name = "Class of call"),
    shares = TRUE,
    title = "How CDD classifies the calls",
    subtitle = paste(
      "Specific: above the model's own threshold. Non-specific: below it.",
      "Superfamily: assigned only to the group of related models."
    )
  )
}

page_agreement <- function(run, layout) {
  calls <- mutate(
    run$calls,
    Method_Agreement = factor(Method_Agreement, CALL_AGREEMENT)
  )
  page_stacked(
    calls, layout, "Method_Agreement",
    scale_fill_agreement(name = "Agreement", limits = CALL_AGREEMENT),
    shares = TRUE,
    title = "Do the two methods agree?",
    subtitle = paste(
      "Share of CDD calls that HMMER confirms: the same family on the same hit",
      "at an overlapping position."
    )
  )
}

page_scores <- function(run, layout) {
  calls <- laid_out(run$calls, layout) |>
    mutate(Score = pmin(-log10(Evalue), 60))
  plot <- ggplot(calls, aes(x = Score, y = Species)) +
    geom_boxplot(
      fill = RAMP_LIGHT, colour = RAMP_DARK, linewidth = 0.35, width = 0.6,
      outlier.size = 0.3, outlier.alpha = 0.3
    ) +
    labs(
      title = "How strong the domain calls are",
      subtitle = paste(
        "Each call's E-value as -log10, capped at 60.", "Further right is stronger."
      ),
      x = "-log10(E-value)"
    )
  species_page(species_by_family(plot), layout)
}

page_queries <- function(run, layout) {
  found <- query_support(run$instances)
  support <- laid_out(found, layout) |>
    mutate(Found_By = factor(N_Queries, levels = rev(seq_along(layout$queries))))
  plot <- ggplot(support, aes(x = N, y = Species, fill = Found_By)) +
    geom_col(width = 0.7, show.legend = TRUE) +
    scale_fill_fixed(
      ordinal_colours(levels(support$Found_By)), name = "Queries that found the domain"
    ) +
    scale_x_continuous(labels = count_labels, expand = expansion(mult = c(0, 0.04))) +
    labs(
      title = "What the queries found together",
      subtitle = paste(
        "Domain instances by the number of query proteins that led to them.",
        "A domain found through several queries is counted once."
      ),
      x = NULL
    )
  species_page(species_by_family(plot), layout)
}

# ---- loci ------------------------------------------------------------------------

page_coverage <- function(run, layout) {
  bins <- laid_out(coverage_bins(run$locus_families), layout)
  plot <- ggplot(bins, aes(x = Coverage, y = Species, fill = N)) +
    geom_tile(colour = PAPER, linewidth = 0.4) +
    scale_fill_amount(name = "Loci") +
    scale_x_continuous(labels = scales::label_percent(), breaks = c(0, 0.5, 1)) +
    labs(
      title = "How much of each domain a locus holds",
      subtitle = paste(
        "Share of the best model covered when a locus's exon pieces are put together.",
        "Loci on the right hold the whole domain."
      ),
      x = "Model coverage"
    ) +
    family_panels() +
    scale_y_discrete(drop = FALSE) +
    labs(y = NULL) +
    theme_rpshunter() +
    italic_rows() +
    theme(panel.grid = element_blank(), legend.key.width = grid::unit(2.2, "cm"))
  species_page(plot, layout)
}

page_copies <- function(run, layout) {
  bins <- copy_bins(run$locus_families, layout$repeat_families) |>
    mutate(
      Species = species_rows(Species, layout),
      Family = factor(Family, layout$repeat_families)
    )
  plot <- ggplot(bins, aes(x = Places, y = Species, fill = N)) +
    geom_tile(colour = PAPER, linewidth = 0.6) +
    tile_numbers(size = 2.8) +
    scale_fill_amount(guide = "none") +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = "Size of repeat arrays",
      subtitle = paste(
        "Loci by the number of separate places a repeat family occurs at.",
        "Families marked `repeat: true` in the config."
      ),
      x = "Places in the locus"
    ) +
    family_panels() +
    scale_y_discrete(drop = FALSE) +
    labs(y = NULL) +
    theme_rpshunter() +
    italic_rows() +
    theme(panel.grid = element_blank())
  species_page(plot, layout, panels = 2 * length(layout$repeat_families))
}

#' The loci drawn on the architecture page, richest first, with a row label
#' that is unique across species.
architecture_rows <- function(run, layout) {
  chosen <- richest_loci(run$loci)
  chosen$Row <- paste(chosen$Species, chosen$Locus)
  chosen$Row <- factor(chosen$Row, levels = rev(chosen$Row))
  chosen$Species_Name <- factor(
    layout$species_names[match(chosen$Species, layout$species)], layout$species_names
  )
  chosen$Length_Kb <- chosen$Length / 1000
  chosen
}

page_architecture <- function(run, layout) {
  title <- "Architecture of the loci carrying most families"
  rows <- architecture_rows(run, layout)
  if (nrow(rows) == 0) {
    return(empty_page(title, "No locus carries more than one family."))
  }
  tracks <- locus_tracks(rows, run$instances) |>
    inner_join(
      select(rows, Species, Locus, Row, Species_Name),
      by = c("Species", "Locus")
    ) |>
    mutate(Family = family_columns(Family, layout))
  plot <- ggplot(rows, aes(y = Row)) +
    geom_segment(
      aes(x = 0, xend = Length_Kb, yend = Row),
      colour = GREY, linewidth = 0.4
    ) +
    geom_segment(
      data = tracks, aes(x = From, xend = To, yend = Row, colour = Family),
      linewidth = 4, show.legend = TRUE
    ) +
    scale_colour_manual(
      values = categorical_colours(layout$families), labels = words,
      limits = layout$families, drop = FALSE, name = "Domain family"
    ) +
    scale_y_discrete(labels = function(row) sub("^\\S+ ", "", row)) +
    facet_wrap(vars(Species_Name), ncol = 1, scales = "free") +
    labs(
      title = title,
      subtitle = paste(
        "Domains at their own position along each locus, exons apart.",
        "Up to eight loci per species, those with most families first."
      ),
      x = "Position in the locus (kb)", y = NULL
    ) +
    theme_rpshunter() +
    theme(
      strip.text = element_text(face = "bold.italic"),
      panel.grid.major.y = element_blank()
    )
  panels <- n_distinct(rows$Species)
  page(plot, page_size(rows = nrow(rows) + 5 * panels, panels = 3))
}
