---
title: Usage & Deployment
author: Tom Stanton
comments: true
tags: [tutorials, web, deployment]
icon: lucide/server
categories:
  - Kaptive-Web
---

Kaptive-Web is designed to be completely flexible. Whether you are a solo researcher running it locally, or a sysadmin deploying it on an HPC or Cloud cluster, there is an execution method for you.

## 💻 Local Execution

### 1. The Native CLI Hook
The easiest way to run the local development server is using the native CLI command installed by the package. This starts the Uvicorn server on `127.0.0.1:8000`.

```bash
kaptive-web
```

If you are using `uv` to manage environments:
```bash
uv run kaptive-web
```

### 2. The Task Runner (`just`)
If you have [just](https://github.com/casey/just) installed, the project comes with a `justfile` containing handy aliases for common tasks. Simply type `just` to see all commands!

- `just serve` - Syncs dependencies and runs the development server.
- `just test` - Runs the pytest suite.
- `just clean` - Cleans up pycache and virtual environments.

## 🐳 Container Deployments

For production deployments or HPC environments, containerization is fully supported.

### Docker & Docker Compose
The `justfile` provides commands for Docker-based deployment:

- `just docker-build`: Builds the standalone `kaptive-web:latest` image.
- `just docker-serve`: Spins up Kaptive-Web and a Caddy reverse proxy using Docker Compose in detached mode.
- `just docker-stop`: Tears down the compose stack.

### Singularity & Apptainer
For HPC environments where Docker is not permitted, you can build Singularity or Apptainer images (`.sif` files).

- `just singularity-build`: Builds via `Singularity.def`
- `just apptainer-build`: Builds via `Apptainer.def`

Once built, you can run the container directly:
```bash
apptainer run kaptive-web.sif
```
