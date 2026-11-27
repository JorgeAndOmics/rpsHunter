# Locate the scripts from the tests, wherever testthat is started from.
SCRIPTS <- normalizePath(file.path("..", "..", "scripts"))
source_script <- function(...) source(file.path(SCRIPTS, ...), local = FALSE)
