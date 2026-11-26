# lintr configuration for the R sources (read by `make lint-r`).
#
# lintr's defaults, with these changes:
# - line length 88, the same limit ruff uses for Python (pyproject.toml);
# - cyclomatic complexity at most 8 per function, the house limit that qlty
#   applies to Python (.qlty/qlty.toml), which cannot measure R;
# - object_usage_linter off: the scripts name data frame columns bare (dplyr,
#   ggplot2), so it reports "no visible binding" for every column;
# - names: snake_case, plus SNAKE_CASE for file-level constants.
#
# LEGACY LIST. The scripts under `exclusions` predate the gate and are being
# rewritten stage by stage, tests first. A whole file is left out until it is
# rewritten, instead of silencing its findings line by line. The list only
# shrinks: delete an entry in the commit that rewrites the file, never add one.
linters <- linters_with_defaults(
  line_length_linter = line_length_linter(88L),
  cyclocomp_linter = cyclocomp_linter(complexity_limit = 8L),
  object_usage_linter = NULL,
  object_name_linter = object_name_linter(styles = c("snake_case", "SNAKE_CASE", "symbols"))
)
exclusions <- list(
  "workflow/scripts/orf_analyser.R"
)
