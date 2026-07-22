# Kaptive Project Justfile
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
test: install
    uv run pytest tests/

# Format Python code
fmt:
    uvx run black .

# Check Python formatting
fmt-check:
    uvx run black --check .

# Lint Python code
lint:
    uvx run ruff check .
    uvx run ty check .

# Run the full CI pipeline locally (format check, lint, test)
ci: fmt-check lint test

# Build and serve the documentation locally
docs-build: install
    uv run scripts/generate_cli_docs.py
    uv run --group docs zensical build

# Test if documentation can be built without warnings or errors
docs-test: install
    uv run scripts/generate_cli_docs.py
    uv run --group docs zensical build -s

serve: install
    uv run scripts/generate_cli_docs.py
    uv run --group docs zensical serve

# Build the Python package
build:
    uv build

# Publish the Python package to PyPI
publish:
    uv publish
