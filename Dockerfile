FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

# Set working directory
WORKDIR /app

# Enable Python bytecode compilation for faster startup
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Copy the dependency manifests
COPY pyproject.toml uv.lock ./

# Install dependencies (without the project itself) to cache the layer
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --extra plot --extra json --no-install-project

# Copy the rest of the project source
COPY src ./src
COPY README.md LICENSE ./

# Install the project
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --extra plot --extra json

# Add the virtual environment to the PATH so the 'kaptive' command is available globally
ENV PATH="/app/.venv/bin:$PATH"

# Set the default entrypoint
ENTRYPOINT ["kaptive"]
CMD ["--help"]
