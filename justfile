# Kaptive Project Justfile
# Run `just` to see all available commands

set shell := ["bash", "-c"]

# Show available commands
default:
    @just --list

# Sync python dependencies and create the virtual environment using `uv`
sync:
    uv sync

# Run Python pytest suite
test: sync
    uv run pytest tests/

# Clean Python virtual environments
clean:
    rm -rf .venv
    find . -type d -name "__pycache__" -exec rm -rf {} +
    find . -type d -name ".pytest_cache" -exec rm -rf {} +

# Boot the Marimo interactive benchmarking server
notebook:
    uv run --group dev marimo edit tests/notebook.py