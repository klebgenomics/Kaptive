# Variables
IDX_FILE := docs/index.md
CLI_FILE := docs/cli.md
LIBNAME  := eris

# Phony targets aren't real files
.PHONY: all docs clean build dev sync test

# Default target when you just run `make`
all: docs

build:
	@echo "🏗️ Building $(LIBNAME) package for distribution using uv..."
	uv build
	@echo "✅ Build complete! Artifacts are in the dist/ directory."

dev:
	@echo "🚧 Installing $(LIBNAME) and dev dependencies..."
	uv sync --all-groups
	@echo "✅ Development installation complete!"

sync:
	@echo "💫 Syncing environment dependencies..."
	uv sync
	@echo "✅ Environment synced!"

#test:
#	@echo "🧪 Running pytest..."
#	uv run --group test pytest -v
#	@echo "✅ Tests complete!"

# The master docs build target
docs: $(CLI_FILE) | $(IDX_FILE)
	@echo "🚀 Building static site..."
	uv run --group docs zensical build --clean
	@echo "✅ Documentation built successfully!"

# Generate the CLI markdown
$(CLI_FILE): $(IDX_FILE)
	@echo "💻 Generating CLI Markdown from $(LIBNAME)..."
	@echo "# CLI Reference" > $@
	@echo "" >> $@
	@echo "$(LIBNAME) provides a highly optimized command-line interface for running the graph traversal pipeline." >> $@
	@echo "" >> $@
	@echo '```shell' >> $@
	uv run $(LIBNAME) -h >> $@
	@echo '```' >> $@
	@echo "✅ Successfully generated CLI documentation at $@"

$(IDX_FILE): README.md
	@echo "📂 Copying README to index..."
	cp README.md $@
	@echo "✅ Completed successfully!"

# Clean target to wipe generated docs and Python build artifacts
clean:
	@echo "🧹 Cleaning up generated documentation and build artifacts..."
	rm -rf $(CLI_FILE) dist/ build/ *.egg-info .pytest_cache
	@echo "✅ Clean complete."