source_script("utils", "log.R")

test_that("a line has time, level, step, scope and message", {
  log_job("figures", "all")
  expect_equal(
    log_line("OK", "3 PDF files", time = "14:02:40"),
    "14:02:40 OK figures all | 3 PDF files"
  )
})

test_that("a message with line breaks stays one line", {
  log_job("figures")
  expect_equal(
    log_line("WARN", "first\n  second", time = "09:00:00"),
    "09:00:00 WARN figures all | first second"
  )
})

test_that("the line shape is the one log.py reads back", {
  log_job("figures", "query_1/Homo_sapiens")
  line <- log_line("INFO", "x | y", time = "01:02:03")
  shape <- "^\\d\\d:\\d\\d:\\d\\d (DEBUG|INFO|OK|WARN|ERROR) \\S+ \\S+ \\| .*$"
  expect_match(line, shape, perl = TRUE)
})

test_that("counts carry thousands separators", {
  expect_equal(log_count(33721), "33,721")
  expect_equal(log_count(7), "7")
})
