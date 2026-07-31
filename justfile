# Kaptive Project Justfile
# Run `just` to see all available commands

set shell := ["bash", "-uc"]

# Show available commands
default:
    @just --list

# Clean caches and build artifacts
clean:
    rm -rf site .cache .ruff_cache docs/api build dist
    find . -type d -name "__pycache__" -exec rm -rf {} +
    find . -type d -name ".pytest_cache" -exec rm -rf {} +
    find . -type f -name "*.pyc" -delete

# Deep clean including the virtual environment
clean-all: clean
    rm -rf .venv

# Sync the virtual environment
install:
    uv sync

# Update dependencies and lockfile
update:
    uv lock --upgrade
    uv sync

# Run the Kaptive CLI (e.g. just run db list)
run *args:
    uv run kaptive {{args}}

# Run the test suite
test:
    uv run pytest tests/

# Format all Python code
fmt:
    uvx ruff format .

# Check if code is formatted without modifying files
fmt-check:
    uvx ruff format --check .

# Lint Python code and auto-fix safe errors
lint:
    uvx ruff check --fix .

# Static type-check Python code
type-check:
    uvx ty check .

# Run all quality checks at once (ideal for local pre-commit testing)
check-all: fmt-check lint type-check

# Run the full CI pipeline locally
ci: check-all test

# Prepare documentation by generating markdown files
prep-docs: install
    uv run scripts/generate_cli_docs.py
    uv run scripts/generate_release_notes.py
    uv run scripts/fetch_announcement.py
    uv run scripts/generate_api_docs.py

# Build the documentation locally
docs-build: prep-docs
    uv run --group docs zensical build

# Translate documentation locally (Requires DEEPL_API_KEY environment variable)
translate: prep-docs
    uv run scripts/translate_docs.py

# Test if documentation can be built without warnings or errors
docs-test: prep-docs
    uv run --group docs zensical build -s

# Serve the documentation locally (English only, with live reload)
serve: prep-docs
    uv run --group docs zensical serve

# Serve the fully built multi-language site locally (no live reload)
serve-all: docs-build
    mkdir -p .local_serve
    ln -sfn ../site .local_serve/Kaptive
    echo "Serving at http://localhost:8000/Kaptive/"
    python3 -m http.server -d .local_serve 8000

# Build the Python package
build: clean
    uv build

# Publish the Python package to PyPI
publish: build
    uv publish
