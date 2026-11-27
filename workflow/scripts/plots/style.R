# The visual system of every rpsHunter figure: colour, words, type, page size
# and output. Plotting code asks this file instead of deciding for itself, and
# never writes a colour value of its own (a test checks). The reasons are in
# the visual style page of the documentation.

suppressPackageStartupMessages(library(ggplot2))

# ---- colour ------------------------------------------------------------------

# Okabe and Ito's colour-blind-safe palette (Color Universal Design, 2008).
PALETTE <- c(
  orange = "#E69F00", sky = "#56B4E9", green = "#009E73", yellow = "#F0E442",
  blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7", black = "#000000"
)
INK <- "#222222"
PAPER <- "#FFFFFF"
SOFT_INK <- "#5A5A5A"
GREY <- "#BBBBBB"
LIGHT_GREY <- "#E6E6E6"
RAMP_LIGHT <- "#EAF3F9"
RAMP_DARK <- PALETTE[["blue"]]

# Fixed meanings: one colour, one meaning, on every page.
EVIDENCE_LEVELS <- c("standard", "borderline", "rescue")
HIT_TYPE_LEVELS <- c("Specific", "Non-specific", "Superfamily")
AGREEMENT_COLOURS <- c(
  both = PALETTE[["green"]], cdd_only = PALETTE[["orange"]],
  hmmer_only = PALETTE[["sky"]], not_compared = GREY
)

#' Colours for an ordered set, strongest first: a dark-to-light ramp.
ordinal_colours <- function(levels) {
  steps <- grDevices::colorRampPalette(c(RAMP_DARK, RAMP_LIGHT))(length(levels) + 1L)
  stats::setNames(steps[seq_along(levels)], levels)
}

#' Colours for an open-ended set (families, queries), in the order given.
#' Past the palette's hues, later levels take lighter tints of the same hues.
categorical_colours <- function(levels) {
  hues <- unname(PALETTE[c(
    "blue", "orange", "green", "purple", "vermillion", "sky", "yellow"
  )])
  round_of <- (seq_along(levels) - 1L) %/% length(hues)
  base <- hues[(seq_along(levels) - 1L) %% length(hues) + 1L]
  tinted <- mapply(tint, base, round_of / (max(round_of) + 1), USE.NAMES = FALSE)
  stats::setNames(tinted, levels)
}

#' Mix a colour with white: 0 leaves it, 1 gives white.
tint <- function(colour, share) {
  rgb <- grDevices::col2rgb(colour)[, 1]
  mixed <- round(rgb + (255 - rgb) * share)
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

#' A fill for a fixed set of levels: every level keeps its key in the legend,
#' present on the page or not.
scale_fill_fixed <- function(values, labels = waiver(), limits = names(values), ...) {
  scale_fill_manual(
    values = values, labels = labels, limits = limits, drop = FALSE, ...
  )
}

scale_fill_evidence <- function(...) {
  scale_fill_fixed(ordinal_colours(EVIDENCE_LEVELS), labels = words, ...)
}

scale_fill_hit_type <- function(...) {
  scale_fill_fixed(ordinal_colours(HIT_TYPE_LEVELS), ...)
}

scale_fill_agreement <- function(...) {
  scale_fill_fixed(AGREEMENT_COLOURS, labels = words, ...)
}

#' A continuous fill for counts, light to dark, on a scale that keeps small
#' counts visible beside large ones.
scale_fill_amount <- function(...) {
  scale_fill_gradient(
    low = RAMP_LIGHT, high = RAMP_DARK, na.value = LIGHT_GREY, trans = "log1p",
    breaks = c(1, 10, 100, 1000, 10000), labels = count_labels, ...
  )
}

#' Counts printed on tiles filled by `scale_fill_amount`: light text on the
#' dark end, dark text elsewhere. `N` is the column holding the count.
tile_numbers <- function(size = 3) {
  list(
    geom_text(
      aes(label = count_labels(N), colour = log1p(N) > 0.75 * log1p(max(N, 1))),
      family = FONT, size = size, show.legend = FALSE
    ),
    scale_colour_manual(
      values = c(`TRUE` = PAPER, `FALSE` = INK), limits = c("TRUE", "FALSE")
    )
  )
}

#' Panel titles as words.
word_labeller <- function() as_labeller(words)

#' A page that says why there is nothing to draw.
empty_page <- function(title, message) {
  plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = message, family = FONT, colour = SOFT_INK) +
    labs(title = title) +
    theme_rpshunter() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          panel.grid = element_blank())
  page(plot, c(width = 9, height = 4))
}

# ---- words -------------------------------------------------------------------

WORDS <- c(
  standard = "Standard (E ≤ 0.01)", borderline = "Borderline (E ≤ 1)",
  rescue = "Rescue (E > 1)", both = "Both methods", cdd_only = "CDD only",
  hmmer_only = "HMMER only", not_compared = "Not compared"
)

#' Plain words for identifiers; anything unlisted has its underscores as spaces.
words <- function(identifiers) {
  known <- WORDS[as.character(identifiers)]
  fallback <- gsub("_", " ", as.character(identifiers), fixed = TRUE)
  unname(ifelse(is.na(known), fallback, known))
}

#' Whole-number counts with thousands separators.
count_labels <- function(values) format(values, big.mark = ",", scientific = FALSE)

# ---- type and theme --------------------------------------------------------------

FONT <- "IBM Plex Sans"

#' The one theme. Species names on an axis are set in italics by the pages.
theme_rpshunter <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = FONT) +
    theme(
      text = element_text(colour = INK),
      plot.title = element_text(face = "bold", size = rel(1.3), hjust = 0),
      plot.subtitle = element_text(colour = SOFT_INK, hjust = 0, lineheight = 1.15),
      plot.title.position = "plot",
      plot.caption = element_text(colour = SOFT_INK, hjust = 0),
      plot.caption.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = LIGHT_GREY, linewidth = 0.3),
      strip.text = element_text(face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.justification = "left",
      legend.title = element_text(face = "bold"),
      panel.spacing.x = grid::unit(1.6, "lines"),
      plot.margin = margin(14, 16, 12, 14)
    )
}

#' Italic row labels, for axes that list species.
italic_rows <- function() theme(axis.text.y = element_text(face = "italic"))

# ---- pages ---------------------------------------------------------------------

MAX_INCHES <- 60

#' Page size in inches for a figure with `rows` rows and `panels` panels across.
#' A page grows with what it holds, so 100 species are as readable as 5, up to
#' a size every viewer opens.
page_size <- function(rows, panels = 1L, row_height = 0.24, panel_width = 2.2) {
  height <- 2.6 + row_height * max(rows, 1L)
  width <- 3.4 + panel_width * max(panels, 1L)
  c(width = min(max(width, 7), MAX_INCHES), height = min(max(height, 4), MAX_INCHES))
}

#' A page: a plot with the size it should be drawn at.
page <- function(plot, size) list(plot = plot, size = size)

#' The key page of a PDF: what each fixed colour means.
key_page <- function(title, families) {
  keys <- rbind(
    key_rows("Evidence tier of a CDD call", ordinal_colours(EVIDENCE_LEVELS)),
    key_rows("How CDD classifies a call", ordinal_colours(HIT_TYPE_LEVELS)),
    key_rows("Method agreement", AGREEMENT_COLOURS),
    key_rows("Domain family (config order)", categorical_colours(families))
  )
  keys$group <- factor(keys$group, levels = unique(keys$group))
  plot <- ggplot(keys, aes(x = 0, y = stats::reorder(label, -order), fill = colour)) +
    geom_tile(width = 0.5, height = 0.7) +
    scale_fill_identity() +
    facet_wrap(vars(group), scales = "free_y", ncol = 2) +
    scale_x_continuous(limits = c(-0.3, 0.3)) +
    labs(
      title = title,
      subtitle = paste(
        "Colours keep one meaning on every page.",
        "Darker means stronger within an ordered set."
      ),
      x = NULL, y = NULL
    ) +
    theme_rpshunter() +
    theme(axis.text.x = element_blank(), panel.grid = element_blank())
  page(plot, c(width = 9, height = 5.5))
}

#' Rows of the key: one per level of a colour map.
key_rows <- function(group, colours) {
  data.frame(
    group = group, label = words(names(colours)), colour = unname(colours),
    order = seq_along(colours)
  )
}

#' Write pages to one vector PDF with the font embedded.
#' A PDF device has one sheet size: the sheet takes the largest page, and every
#' page fills it. The pages of one file are sized alike (same species rows), so
#' none is stretched far.
write_pdf <- function(pages, path) {
  sizes <- vapply(pages, function(each) each$size, numeric(2))
  grDevices::cairo_pdf(
    path, width = max(sizes["width", ]), height = max(sizes["height", ]),
    onefile = TRUE
  )
  on.exit(grDevices::dev.off())
  for (each in pages) print(each$plot)
  invisible(path)
}
