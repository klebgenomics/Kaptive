# Kaptive Project Justfile
# Run `just` to see all available commands

set shell := ["bash", "-c"]

# Show available commands
default:
    @just --list

# Sync python dependencies and create the virtual environment using `uv`
sync:
    uv sync

# Generate CLI documentation from the binary help output
doc-cli: build
    mkdir -p docs
    uv run python docs/generate_cli_docs.py

# Generate benchmark documentation
doc-bench: build
    mkdir -p docs
    echo "# Benchmarks" > docs/benchmarks.md
    echo "\`\`\`text" >> docs/benchmarks.md
    cargo bench --workspace >> docs/benchmarks.md || true
    echo "\`\`\`" >> docs/benchmarks.md

# Build all documentation and generate Zensical site
docs: build doc-cli doc-bench
    cargo doc --no-deps --workspace
    mkdir -p docs/rust-api
    cp -r target/doc/* docs/rust-api/
    cp README.md docs/index.md
    uv run zensical build

# Serve the documentation site locally
docs-serve: docs
    uv run zensical serve

# Deploy the documentation to GitHub Pages
docs-deploy: docs
    uv run zensical gh-deploy --force

# Run Python pytest suite
test: sync
    uv run pytest tests/

# Clean Python virtual environments
clean:
    rm -rf .venv
    find . -type d -name "__pycache__" -exec rm -rf {} +
    find . -type d -name ".pytest_cache" -exec rm -rf {} +