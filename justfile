# pyfgs Project Justfile
# Run `just` to see all available commands

set shell := ["bash", "-uc"]

# Show available commands
default:
    @just --list

# Clean build artifacts and cache
clean:
    rm -rf site target/wheels
    # rm -rf .venv
    find . -type d -name "__pycache__" -exec rm -rf {} +
    find . -type d -name ".pytest_cache" -exec rm -rf {} +
    find . -type d -name ".karva_cache" -exec rm -rf {} +

# Install dependencies, compile Rust extension in dev mode, and sync stubs
install: clean
    uv sync --all-groups
    cargo run --bin stub_gen --no-default-features
    uv run python scripts/sync_stubs.py
    uvx ruff format python/pyfgs/__init__.py python/pyfgs/__init__.pyi python/pyfgs/_pyfgs/__init__.pyi

# Run the test suite using karva
test +args="": install
    uvx --with . --with pytest karva test tests/ {{args}}

# Run tests with coverage using karva
test-cov +args="": install
    uvx --with . --with pytest karva test --cov tests/ {{args}}

# Run benchmarks
benchmark: install
    uv run --group benches pytest benchmarks/bench_accuracy.py --benchmark-json=benchmark-data.json

# --- Formatting & Linting ---

# Format both Python and Rust code
fmt: fmt-python fmt-rust

# Check Python and Rust formatting
fmt-check: fmt-check-python fmt-check-rust

# Lint both Python and Rust code
lint: lint-python lint-rust

# Format Python code
fmt-python:
    uvx ruff format .

# Check Python formatting
fmt-check-python:
    uvx ruff format --check .

# Lint Python code
lint-python:
    -uvx ruff check .
    -uv run ty check .

# Format Rust code
fmt-rust:
    cargo fmt

# Check Rust formatting
fmt-check-rust:
    cargo fmt --all -- --check

# Lint Rust code
lint-rust:
    cargo clippy --workspace --all-targets --all-features -- -D warnings

# --- CI & Build ---

# Run the full CI pipeline locally (format check, lint, test)
ci: install fmt-check lint test-cov

# Build release wheels using maturin
build-wheels:
    rm -rf target/wheels
    uvx maturin build --release

# Set version in Cargo.toml (useful in CI or releases)
set-version VERSION:
    uv run python -c "import re; content = open('Cargo.toml').read(); content = re.sub(r'^version\s*=\s*\".*\"', 'version = \"' + '{{VERSION}}'.lstrip('v') + '\"', content, flags=re.MULTILINE); open('Cargo.toml', 'w').write(content)"

# Generate release notes from GitHub API
generate-releases:
    uv run python scripts/generate_release_notes.py

# --- Documentation ---

prep-docs: install
    uv run scripts/generate_cli_docs.py
    uv run python scripts/generate_api_docs.py
    uv run python scripts/fetch_announcement.py

# Build the documentation locally
docs-build: prep-docs
    uv run --group docs zensical build -c

# Test if documentation can be built without warnings or errors
docs-test: prep-docs
    uv run --group docs zensical build -c -s

serve: prep-docs
    uv run --group docs zensical serve

# Build the Python package
build:
    uv build

# Publish the Python package to PyPI
publish:
    uv publish
