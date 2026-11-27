source_script("plots", "tables.R")

layout <- list(
  species = c("B_two", "A_one", "C_none"),
  species_names = c("B two", "A one", "C none"),
  families = c("first", "second"), repeat_families = "second",
  queries = c("NP_1.1", "NP_2.1")
)
run <- list(
  hit_species = c(rep("A_one", 5), rep("B_two", 3)),
  calls = data.frame(
    Species = c("A_one", "A_one", "A_one", "B_two"),
    Tag = c("t1", "t1", "t2", "t3"),
    Family = c("first", "second", "first", "first")
  ),
  loci = data.frame(
    Species = c("A_one", "A_one", "B_two"),
    Locus = c("s:1-9", "s:50-90", "s:1-5"),
    Sequence = "s", Start = c(1, 50, 1), End = c(9000, 90000, 5000),
    N_Hits = c(3L, 1L, 2L), N_Families = c(2L, 0L, 1L),
    Best_Bitscore = c(90, 20, 50)
  ),
  locus_families = data.frame(
    Species = c("A_one", "A_one", "B_two"),
    Locus = c("s:1-9", "s:1-9", "s:1-5"),
    Family = c("first", "second", "first"),
    Coverage = c(1, 0.34, 0.05), Copies = c(1L, 7L, 1L)
  )
)

test_that("tip labels are read from a Newick tree in written order", {
  tree <- "((Homo_sapiens:1,'Mus musculus':2):0.5,Desmodus_rotundus:3);"
  tips <- c("Homo_sapiens", "Mus musculus", "Desmodus_rotundus")
  expect_equal(newick_tips(tree), tips)
})

test_that("species follow the tree, by key or display name, the rest last", {
  ordered <- tree_order(
    keys = c("Mus_musculus", "Extra", "Homo_sapiens"),
    names = c("Mus musculus", "An extra one", "Homo sapiens"),
    tips = c("Homo sapiens", "mus_musculus")
  )
  expect_equal(ordered, c("Homo_sapiens", "Mus_musculus", "Extra"))
})

test_that("the layout comes from the config, in config order without a tree", {
  config <- list(
    species = list(B_two = list(name = "B two"), A_one = list(name = "A one")),
    domains = list(first = list(cdd = "x"), second = list(cdd = "y", `repeat` = TRUE)),
    queries = list(NP_1.1 = "q")
  )
  found <- run_layout(config)
  expect_equal(found$species, c("B_two", "A_one"))
  expect_equal(found$species_names, c("B two", "A one"))
  expect_equal(found$families, c("first", "second"))
  expect_equal(found$repeat_families, "second")
})

test_that("every configured species keeps its row, the first on top", {
  rows <- species_rows(c("A_one", "B_two"), layout)
  expect_equal(levels(rows), c("C none", "A one", "B two"))
})

test_that("the funnel counts each step for every species, zero included", {
  counts <- funnel_counts(run, layout)
  of <- function(species, step) {
    counts$N[counts$Species == species & counts$Step == step]
  }
  expect_equal(of("A_one", "tBLASTn hits"), 5)
  expect_equal(of("A_one", "Selected hits"), 4)
  expect_equal(of("A_one", "Hits with a domain"), 2)
  expect_equal(of("A_one", "Loci"), 2)
  expect_equal(of("A_one", "Loci with a family"), 1)
  expect_equal(of("C_none", "Loci"), 0)
  expect_equal(nrow(counts), 3 * length(FUNNEL_STEPS))
})

test_that("family presence counts loci per family and loci with every family", {
  presence <- family_presence(run, layout)
  of <- function(species, family) {
    presence$N[presence$Species == species & presence$Family == family]
  }
  expect_equal(of("A_one", "first"), 1)
  expect_equal(of("A_one", ALL_FAMILIES), 1)
  expect_equal(of("B_two", ALL_FAMILIES), 0)
  expect_equal(of("C_none", "second"), 0)
  expect_equal(nrow(presence), 3 * 3)
})

test_that("instances are counted by the number of queries that found them", {
  instances <- data.frame(
    Species = "A_one", Family = "first",
    Queries = c("NP_1.1", "NP_1.1,NP_2.1", "NP_2.1")
  )
  support <- query_support(instances)
  expect_equal(support$N[support$N_Queries == 1], 2)
  expect_equal(support$N[support$N_Queries == 2], 1)
})

test_that("coverage falls into tenths, full coverage in the last one", {
  bins <- coverage_bins(run$locus_families)
  expect_equal(sort(bins$Coverage), c(0.05, 0.35, 0.95))
})

test_that("copies are binned for repeat families only", {
  bins <- copy_bins(run$locus_families, layout$repeat_families)
  expect_equal(as.character(bins$Family), "second")
  expect_equal(as.character(bins$Places), "6 to 10")
})

test_that("the richest loci carry at least two families, best first, capped", {
  chosen <- richest_loci(run$loci, per_species = 8L)
  expect_equal(chosen$Locus, "s:1-9")
})

test_that("domains are placed inside their locus in kilobases from its start", {
  chosen <- richest_loci(run$loci)
  instances <- data.frame(
    Species = "A_one", Sequence = "s", Start = c(1001, 20000), End = c(2000, 20100),
    Family = c("first", "second")
  )
  tracks <- locus_tracks(chosen, instances)
  expect_equal(nrow(tracks), 1) # the second lies outside the locus
  expect_equal(tracks$From, 1)
  expect_equal(tracks$To, 2)
})
