source_script("plots", "style.R")
source_script("plots", "tables.R")
source_script("plots", "pages.R")

# A run with `n` species, built from nothing: enough of every table for every
# page, with one species left without any result.
synthetic_run <- function(n, families = c("first", "second", "third")) {
  species <- sprintf("Genus_species%03d", seq_len(n))
  with_data <- head(species, n - 1)
  layout <- list(
    species = species, species_names = gsub("_", " ", species),
    families = families, repeat_families = "third", queries = c("NP_1.1", "NP_2.1")
  )
  grid <- expand.grid(
    Species = with_data, Family = families, Rank = 1:4, stringsAsFactors = FALSE
  )
  loci <- data.frame(
    Species = rep(with_data, each = 3), Locus = paste0("s:", 1:3),
    Sequence = "s", Start = rep(c(1, 50001, 90001), length(with_data)),
    End = rep(c(20000, 62000, 99000), length(with_data)), Length = 20000,
    N_Hits = 4L, N_Families = rep(c(3L, 2L, 0L), length(with_data)),
    Best_Bitscore = rep(c(90, 50, 20), length(with_data))
  )
  run <- list(
    hit_species = rep(with_data, each = 12),
    calls = transform(
      grid,
      Tag = paste0("t", Rank), Evalue = 10^-(Rank * 3),
      Evidence = ifelse(Rank == 1, "borderline", "standard"),
      Hit_Type = ifelse(Rank < 3, "Specific", "Superfamily"),
      Method_Agreement = ifelse(Rank == 2, "cdd_only", "both")
    ),
    instances = transform(
      grid,
      Sequence = "s", Start = Rank * 3000, End = Rank * 3000 + 300,
      Queries = ifelse(Rank == 4, "NP_1.1,NP_2.1", "NP_1.1")
    ),
    loci = loci,
    locus_families = data.frame(
      Species = rep(with_data, each = 3), Locus = "s:1", Family = families,
      Coverage = rep(c(1, 0.42, 0.08), length(with_data)),
      Copies = rep(c(1L, 3L, 14L), length(with_data))
    )
  )
  list(run = run, layout = layout)
}

all_pages <- function(made) {
  with(made, list(
    page_family_presence(run, layout), page_funnel(run, layout),
    page_evidence(run, layout), page_scores(run, layout),
    page_hit_type(run, layout), page_agreement(run, layout),
    page_queries(run, layout), page_coverage(run, layout),
    page_copies(run, layout), page_architecture(run, layout)
  ))
}

row_labels <- function(plot) {
  built <- ggplot_build(plot)
  built$layout$panel_params[[1]]$y$get_labels()
}

test_that("every page draws at 5, 30 and 100 species", {
  skip_if_not(capabilities("cairo"))
  for (n in c(5, 30, 100)) {
    pages <- all_pages(synthetic_run(n))
    path <- tempfile(fileext = ".pdf")
    expect_no_error(write_pdf(pages, path))
    expect_gt(file.size(path), 5000)
  }
})

test_that("pages with species on rows grow with the number of species", {
  heights <- vapply(c(5, 30, 100), function(n) {
    made <- synthetic_run(n)
    page_funnel(made$run, made$layout)$size[["height"]]
  }, numeric(1))
  expect_true(all(diff(heights) > 0))
  expect_gte(heights[3] / 100, 0.2)
})

test_that("a species with no result keeps its row", {
  made <- synthetic_run(5)
  for (make in list(page_evidence, page_funnel, page_coverage, page_family_presence)) {
    labels <- row_labels(make(made$run, made$layout)$plot)
    expect_true("Genus species005" %in% labels)
    expect_length(labels, 5)
  }
})

test_that("species are listed in layout order, the first on top", {
  made <- synthetic_run(5)
  labels <- row_labels(page_funnel(made$run, made$layout)$plot)
  expect_equal(labels[length(labels)], "Genus species001")
})

test_that("a run without any multi-family locus still gets an architecture page", {
  made <- synthetic_run(5)
  made$run$loci$N_Families <- 1L
  drawn <- page_architecture(made$run, made$layout)
  expect_s3_class(drawn$plot, "ggplot")
})

test_that("a run with no domain at all still draws every page", {
  skip_if_not(capabilities("cairo"))
  made <- synthetic_run(5)
  made$run$calls <- made$run$calls[0, ]
  made$run$instances <- made$run$instances[0, ]
  made$run$locus_families <- made$run$locus_families[0, ]
  made$run$loci$N_Families <- 0L
  expect_no_error(write_pdf(all_pages(made), tempfile(fileext = ".pdf")))
})
