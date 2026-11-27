#!/usr/bin/env Rscript
# Draw the figures of a finished run: four vector PDFs.
#
#   overview.pdf      loci by family, from hits to loci
#   domains.pdf       evidence, strength, CDD's classes, method agreement, queries
#   loci.pdf          model coverage, repeat arrays
#   architecture.pdf  the loci carrying most families, domain by domain
#
# The first three open with a key to the colours. The visual system is in
# plots/style.R, the counts in plots/tables.R and the pages in plots/pages.R.
#
# Usage (called by the Snakefile):
#   Rscript figures.R --results RESULTS --config CONFIG.yaml --out-dir DIR \
#       [--tree NEWICK]
#
# A run with no domain call still gets its PDFs: pages show empty rows.

script_dir <- function() {
  flag <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", flag[1])))
}

load_parts <- function(here) {
  for (part in c("utils/log.R", "plots/style.R", "plots/tables.R", "plots/pages.R")) {
    source(file.path(here, part), local = FALSE)
  }
}

arguments <- function() {
  parser <- argparse::ArgumentParser(description = "Draw the figures of a run.")
  parser$add_argument("--results", required = TRUE)
  parser$add_argument("--config", required = TRUE)
  parser$add_argument("--out-dir", required = TRUE)
  parser$add_argument("--tree", default = NULL)
  parser$parse_args()
}

#' The pages of each PDF, by file name. Every file is always written, so that
#' a run's outputs do not depend on what was found.
pdf_pages <- function(run, layout) {
  key <- function(title) key_page(title, layout$families)
  loci <- list(key("Loci: key to the colours"), page_coverage(run, layout))
  if (length(layout$repeat_families) > 0) {
    loci <- c(loci, list(page_copies(run, layout)))
  }
  list(
    overview.pdf = list(
      key("Overview: key to the colours"),
      page_family_presence(run, layout), page_funnel(run, layout)
    ),
    domains.pdf = list(
      key("Domains: key to the colours"),
      page_evidence(run, layout), page_scores(run, layout),
      page_hit_type(run, layout), page_agreement(run, layout),
      page_queries(run, layout)
    ),
    loci.pdf = loci,
    architecture.pdf = list(page_architecture(run, layout))
  )
}

main <- function() {
  args <- arguments()
  log_job("figures")
  layout <- run_layout(yaml::read_yaml(args$config), args$tree)
  run <- read_run(args$results)
  files <- pdf_pages(run, layout)
  for (name in names(files)) {
    write_pdf(files[[name]], file.path(args$out_dir, name))
  }
  pages <- sum(lengths(files))
  log_ok(sprintf("%d PDF files, %d pages", length(files), pages))
}

if (sys.nframe() == 0) {
  load_parts(script_dir())
  run_main(main)
}
