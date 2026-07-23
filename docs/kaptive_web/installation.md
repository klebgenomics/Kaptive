---
title: Installation
author: Tom Stanton
comments: true
tags: [tutorials, installation, web]
icon: lucide/download
categories:
  - Kaptive-Web
---

You can install `kaptive-web` using standard Python package managers. We highly recommend using `uv` for lightning-fast dependency resolution.

=== "uv (Recommended)"

    ```bash
    # Install globally or in a virtual environment
    uv pip install kaptive-web
    ```

=== "pip"

    ```bash
    pip install kaptive-web
    ```

## :lucide-blocks: Dependencies

The backend is built in modern Python (3.11+) using **FastAPI** and relies on the following core dependencies:

 - [kaptive](https://klebgenomics.github.io/Kaptive/): The core serotyping engine (installed with `plot` and `json` extras).
 - [fastapi](https://fastapi.tiangolo.com/): The web framework.
 - [orjson](https://github.com/ijl/orjson): For blazing fast serialization of NumPy arrays.
 - [aiosqlite](https://aiosqlite.omnilib.dev/en/stable/): For asynchronous SQLite database management.
