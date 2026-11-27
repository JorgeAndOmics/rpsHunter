source_script("plots", "style.R")

# Smallest CIELAB distance between any two colours as seen with each kind of
# colour blindness: the measure of whether a palette stays distinguishable.
min_distance_colour_blind <- function(colours) {
  seen <- list(
    colorspace::deutan(colours), colorspace::protan(colours),
    colorspace::tritan(colours)
  )
  smallest <- vapply(seen, function(set) {
    lab <- farver::decode_colour(set, to = "lab")
    distances <- as.matrix(stats::dist(lab))
    min(distances[upper.tri(distances)])
  }, numeric(1))
  min(smallest)
}

test_that("the categorical hues stay apart under every colour blindness", {
  hues <- unname(categorical_colours(paste0("family_", 1:7)))
  expect_gt(min_distance_colour_blind(hues), 15)
})

test_that("the method-agreement colours stay apart under every colour blindness", {
  expect_gt(min_distance_colour_blind(unname(AGREEMENT_COLOURS)), 15)
})

test_that("an ordered set runs from dark to light", {
  colours <- ordinal_colours(EVIDENCE_LEVELS)
  lightness <- farver::decode_colour(unname(colours), to = "lab")[, "l"]
  expect_named(colours, EVIDENCE_LEVELS)
  expect_true(all(diff(lightness) > 0))
})

test_that("families keep the config order and never share a colour", {
  for (count in c(1, 4, 7, 8, 20)) {
    families <- paste0("family_", seq_len(count))
    colours <- categorical_colours(families)
    expect_named(colours, families)
    expect_equal(anyDuplicated(unname(colours)), 0L)
  }
})

test_that("the first families get the same colour however many follow", {
  few <- categorical_colours(c("a", "b", "c"))
  expect_equal(categorical_colours(c("a", "b", "c", "d"))[1:3], few)
})

test_that("identifiers are shown as words", {
  expect_equal(words("cdd_only"), "CDD only")
  expect_equal(words("zinc_finger"), "zinc finger")
  expect_equal(words(c("both", "SET")), c("Both methods", "SET"))
})

test_that("a page grows with its rows and stays within what viewers open", {
  small <- page_size(rows = 5)
  medium <- page_size(rows = 30)
  large <- page_size(rows = 100)
  expect_lt(small[["height"]], medium[["height"]])
  expect_lt(medium[["height"]], large[["height"]])
  expect_gte(large[["height"]] / 100, 0.2) # a row keeps a readable height
  expect_lte(page_size(rows = 10000, panels = 500)[["height"]], MAX_INCHES)
  expect_lte(page_size(rows = 10000, panels = 500)[["width"]], MAX_INCHES)
})

test_that("pages of different sizes are written to one PDF with the font embedded", {
  skip_if_not(capabilities("cairo"))
  plot <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) +
    geom_point() +
    theme_rpshunter()
  pages <- list(
    key_page("Key", c("first", "second")),
    page(plot, page_size(rows = 30))
  )
  path <- tempfile(fileext = ".pdf")
  write_pdf(pages, path)
  expect_gt(file.size(path), 2000)
  header <- readBin(path, "raw", 5)
  expect_equal(rawToChar(header), "%PDF-")
  expect_true(any(grepl("IBMPlex", readLines(path, warn = FALSE), useBytes = TRUE)))
})

test_that("no plotting code outside the style module writes a colour value", {
  files <- list.files(file.path(SCRIPTS, "plots"), pattern = "\\.R$", full.names = TRUE)
  others <- files[basename(files) != "style.R"]
  for (file in c(others, file.path(SCRIPTS, "figures.R"))) {
    if (!file.exists(file)) next
    lines <- readLines(file)
    expect_false(
      any(grepl("#[0-9A-Fa-f]{6}\\b", lines)),
      info = paste(basename(file), "contains a colour value")
    )
  }
})
