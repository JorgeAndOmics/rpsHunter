.PHONY: help env env-update lint lint-r lint-sh format format-check typecheck test test-py test-r check quality clean

# Every target expects the rpsHunter conda environment to be active
# (`conda activate rpsHunter`), or to be run as `conda run -n rpsHunter make <target>`.

ENV_FILE := data/config/environment.yaml
ENV_LOCK := data/config/environment.lock-linux-64.txt
ENV_NAME := rpsHunter
PY_SRC   := workflow/scripts
PY_TESTS := tests
R_TESTS  := workflow/tests/testthat
# Legacy shell scripts (see the legacy list in pyproject.toml): left out of
# ShellCheck as whole files until rewritten. This list only shrinks.
SH_LEGACY := workflow/scripts/query_downloader.sh
SH_SRC   := $(filter-out $(SH_LEGACY),$(wildcard workflow/scripts/*.sh))

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

# -- environment ------------------------------------------
env: ## Create the conda env exactly as locked (linux-64)
	conda create --yes --name $(ENV_NAME) --file $(ENV_LOCK)

env-update: ## Re-solve the env from environment.yaml (changes versions: see docs)
	conda env update --file $(ENV_FILE) --name $(ENV_NAME) --prune

# -- linting / formatting / typing -----------------------
lint: ## ruff lint on Python sources and tests
	ruff check $(PY_SRC) $(PY_TESTS)

lint-r: ## lintr on R sources (config: .lintr.R)
	Rscript -e 'l <- lintr::lint_dir("workflow"); print(l); quit(status = length(l) > 0)'

lint-sh: ## ShellCheck on the shell scripts
	@if [ -n "$(SH_SRC)" ]; then shellcheck $(SH_SRC); \
	else echo "lint-sh: no shell scripts outside the legacy list"; fi

format: ## Auto-format Python and apply ruff fixes
	ruff format $(PY_SRC) $(PY_TESTS)
	ruff check --fix $(PY_SRC) $(PY_TESTS)

format-check: ## Verify formatting without modifying
	ruff format --check $(PY_SRC) $(PY_TESTS)

typecheck: ## mypy strict on Python sources outside the legacy list
	mypy

# -- testing ---------------------------------------------
test: test-py test-r ## Run Python and R tests

test-py: ## pytest, with the coverage floor from pyproject.toml
	pytest --cov --cov-report=term-missing:skip-covered

# No R test exists yet: the R scripts are all on the legacy list and run top to
# bottom on import. The first testthat files arrive with the first rewritten R
# module; until then this target says so instead of passing silently.
test-r: ## R testthat suite
	@if ls $(R_TESTS)/test-*.R >/dev/null 2>&1; then \
		Rscript -e 'testthat::test_dir("$(R_TESTS)", stop_on_failure = TRUE)'; \
	else echo "test-r: NO R TESTS YET ($(R_TESTS) is empty)"; fi

# -- combined gate ---------------------------------------
check: lint lint-r lint-sh format-check typecheck test ## Run every quality gate (before each commit)

quality: ## Function-size and complexity metrics for the whole tree (qlty)
	qlty metrics --all --functions

# -- cleanup ---------------------------------------------
clean: ## Remove Python and R caches
	rm -rf .pytest_cache .mypy_cache .ruff_cache htmlcov .coverage coverage.xml
	find . -type d -name __pycache__ -not -path './results/*' -exec rm -rf {} + 2>/dev/null || true
	rm -f workflow/scripts/.RData workflow/scripts/.Rhistory
