# One line shape for every message an R job writes, the same as log.py:
#
#   14:02:40 OK figures all | 3 PDF files, 14 pages
#   time     level step    scope  message
#
# Scripts write plain lines; only the launcher adds colour (console.py).
# See docs/console_style.md.

.log_state <- new.env(parent = emptyenv())
.log_state$step <- "r"
.log_state$scope <- "all"

#' Set the step and scope of this job's lines. Call once, first thing in main.
log_job <- function(step, scope = "all") {
  .log_state$step <- step
  .log_state$scope <- scope
  invisible(NULL)
}

#' Format one line. `time` is passed in so that tests are not tied to the clock.
log_line <- function(level, message, time = format(Sys.time(), "%H:%M:%S")) {
  text <- gsub("\\s+", " ", trimws(message))
  sprintf("%s %s %s %s | %s", time, level, .log_state$step, .log_state$scope, text)
}

log_info <- function(message) cat(log_line("INFO", message), "\n", sep = "")
log_warn <- function(message) cat(log_line("WARN", message), "\n", sep = "")
log_ok <- function(message) cat(log_line("OK", message), "\n", sep = "")
log_error <- function(message) cat(log_line("ERROR", message), "\n", sep = "")

#' A count with thousands separators, as messages show them.
log_count <- function(number) format(number, big.mark = ",", scientific = FALSE)

#' Run a script's main; a failure becomes one ERROR line and exit status 1.
run_main <- function(main) {
  status <- tryCatch(
    {
      main()
      0L
    },
    error = function(problem) {
      log_error(conditionMessage(problem))
      1L
    }
  )
  quit(save = "no", status = status)
}
