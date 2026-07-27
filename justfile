# pyfgs Project Justfile
# Run `just` to see all available commands

set shell := ["bash", "-uc"]

# Show available commands
default:
    @just --list

# Clean Python virtual environments
clean:
    rm -rf site
    # rm -rf .venv
    find . -type d -name "__pycache__" -exec rm -rf {} +
    find . -type d -name ".pytest_cache" -exec rm -rf {} +

install: clean
    uv sync

# Run the test suite
test +args="": install
    uv run pytest tests/ {{args}}

# Run tests with coverage
test-cov +args="": install
    uv run --group test pytest -v tests/ --cov=python/pyfgs --cov-report=term-missing --cov-report=xml {{args}}

# Run benchmarks
benchmark: install
    uv run --group benches pytest benchmarks/bench_accuracy.py --benchmark-json=benchmark-data.json

# Format Python code
fmt:
    uv run ruff format .

# Check Python formatting
fmt-check:
    uv run ruff format --check .

# Lint Python code
lint:
    uv run ruff check .
    uv run ty check .

# Run the full CI pipeline locally (format check, lint, test)
ci: fmt-check lint test-cov

cli-docs: install
    uv run scripts/generate_cli_docs.py

prep-docs: cli-docs
    cp README.md docs/index.md
    cp CONTRIBUTING.md docs/contributing.md

# Build the documentation locally
docs-build: prep-docs
    uv run --group docs zensical build

# Test if documentation can be built without warnings or errors
docs-test: prep-docs
    uv run --group docs zensical build -s

serve: prep-docs
    uv run --group docs zensical serve

# Build the Python package
build:
    uv build

# Publish the Python package to PyPI
publish:
    uv publish
