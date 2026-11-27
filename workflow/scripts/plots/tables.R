# What the figures are drawn from: the run's tables, read once, and the counts
# and selections each page needs. No plotting here.

suppressPackageStartupMessages(library(dplyr))

FUNNEL_STEPS <- c(
  "tBLASTn hits", "Selected hits", "Hits with a domain", "Loci", "Loci with a family"
)
COVERAGE_BREAKS <- seq(0, 1, by = 0.1)
COPY_BREAKS <- c(0, 1, 2, 5, 10, 20, Inf)
COPY_LABELS <- c("1", "2", "3 to 5", "6 to 10", "11 to 20", "more than 20")
ALL_FAMILIES <- "Every family"

#' Read the tables of a finished run. Only the columns the pages use are read
#' from the two large tables.
read_run <- function(results) {
  table <- function(folder, name, ...) {
    arrow::read_parquet(file.path(results, "tables", folder, name), ...)
  }
  call_columns <- c(
    "Species", "Tag", "Family", "Evidence", "Evalue", "Hit_Type", "Method_Agreement"
  )
  hits <- table("blast", "aggregate.parquet", col_select = "Species")
  list(
    hit_species = hits$Species,
    calls = table("domains", "calls.parquet", col_select = dplyr::all_of(call_columns)),
    instances = table("domains", "domains.parquet"),
    loci = table("loci", "loci.parquet"),
    locus_families = table("loci", "locus_families.parquet")
  )
}

#' Tip labels of a Newick tree, in the order they are written.
newick_tips <- function(text) {
  labels <- regmatches(text, gregexpr("(?<=[(,])\\s*[^(),:;]+", text, perl = TRUE))[[1]]
  trimws(gsub("'", "", labels, fixed = TRUE))
}

#' Species keys in the order of a tree's tips; species the tree lacks go last.
tree_order <- function(keys, names, tips) {
  plain <- function(labels) tolower(gsub("[_ ]+", " ", labels))
  absent <- length(tips) + 1L
  position <- pmin(
    match(plain(keys), plain(tips), nomatch = absent),
    match(plain(names), plain(tips), nomatch = absent)
  )
  keys[order(position, seq_along(keys))]
}

#' What the config says about how to lay a figure out: species (in tree order
#' if a tree is given, else config order), their display names, the families
#' in config order, and which of them are repeats.
run_layout <- function(config, tree_file = NULL) {
  keys <- names(config$species)
  names <- vapply(config$species, function(entry) entry$name, character(1))
  if (!is.null(tree_file)) {
    tips <- newick_tips(paste(readLines(tree_file, warn = FALSE), collapse = ""))
    keys <- tree_order(keys, names[keys], tips)
  }
  repeats <- vapply(config$domains, function(f) isTRUE(f[["repeat"]]), logical(1))
  list(
    species = keys, species_names = unname(names[keys]),
    families = names(config$domains), repeat_families = names(config$domains)[repeats],
    queries = names(config$queries)
  )
}

#' Species as a factor that keeps every configured species, first one on top.
species_rows <- function(species, layout) {
  factor(species, levels = rev(layout$species), labels = rev(layout$species_names))
}

#' Families as a factor in config order.
family_columns <- function(families, layout, extra = NULL) {
  factor(families, levels = c(layout$families, extra))
}

#' Counts per species at each step from hits to loci with a family.
funnel_counts <- function(run, layout) {
  count_by <- function(species) as.vector(table(factor(species, layout$species)))
  per_step <- list(
    count_by(run$hit_species),
    as.vector(tapply(run$loci$N_Hits, factor(run$loci$Species, layout$species), sum)),
    count_by(distinct(run$calls, Species, Tag)$Species),
    count_by(run$loci$Species),
    count_by(run$loci$Species[run$loci$N_Families > 0])
  )
  data.frame(
    Species = rep(layout$species, times = length(FUNNEL_STEPS)),
    Step = factor(rep(FUNNEL_STEPS, each = length(layout$species)), FUNNEL_STEPS),
    N = replace(unlist(per_step), is.na(unlist(per_step)), 0)
  )
}

#' Loci carrying each family, and loci carrying every family, per species.
family_presence <- function(run, layout) {
  per_family <- count(run$locus_families, Species, Family, name = "N")
  complete <- run$loci |>
    filter(N_Families == length(layout$families)) |>
    count(Species, name = "N") |>
    mutate(Family = ALL_FAMILIES)
  grid <- expand.grid(
    Species = layout$species, Family = c(layout$families, ALL_FAMILIES),
    stringsAsFactors = FALSE
  )
  left_join(grid, bind_rows(per_family, complete), by = c("Species", "Family")) |>
    mutate(N = coalesce(N, 0L))
}

#' Domain instances by how many of the queries found them.
query_support <- function(instances) {
  instances |>
    mutate(N_Queries = lengths(strsplit(Queries, ",", fixed = TRUE))) |>
    count(Species, Family, N_Queries, name = "N")
}

#' Loci per species and family, by tenths of model coverage.
coverage_bins <- function(locus_families) {
  locus_families |>
    mutate(
      Bin = cut(Coverage, COVERAGE_BREAKS, include.lowest = TRUE, labels = FALSE)
    ) |>
    count(Species, Family, Bin, name = "N") |>
    mutate(Coverage = COVERAGE_BREAKS[Bin] + 0.05)
}

#' Loci per species by how many separate places a repeat family occurs at.
copy_bins <- function(locus_families, repeat_families) {
  locus_families |>
    filter(Family %in% repeat_families) |>
    mutate(Places = cut(Copies, COPY_BREAKS, labels = COPY_LABELS)) |>
    count(Species, Family, Places, name = "N")
}

#' The loci carrying most families, at most `per_species` for each species.
richest_loci <- function(loci, per_species = 8L) {
  loci |>
    filter(N_Families >= 2) |>
    group_by(Species) |>
    arrange(desc(N_Families), desc(Best_Bitscore), .by_group = TRUE) |>
    slice_head(n = per_species) |>
    ungroup()
}

#' Domain instances inside the chosen loci, positioned from each locus's start (kb).
locus_tracks <- function(chosen, instances) {
  inside <- inner_join(
    select(chosen, Species, Locus, Sequence, Locus_Start = Start, Locus_End = End),
    select(instances, Species, Sequence, Start, End, Family),
    by = join_by(Species, Sequence, Locus_Start <= Start, Locus_End >= End)
  )
  mutate(
    inside,
    From = (Start - Locus_Start) / 1000, To = (End - Locus_Start + 1) / 1000
  )
}
