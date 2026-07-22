---
title: Installation
author: Tom Stanton
comments: true
tags: [tutorials, installation]
icon: lucide/download
categories:
  - Tutorials
---

=== "pip"

    ```bash
    pip install kaptive
    ```

=== "uv"

    ```bash
    uv add kaptive
    ```

=== "conda"

    ```bash
    conda install bioconda::kaptive
    ```

=== "mamba"

    ```bash
    mamba install bioconda::kaptive
    ```

=== "pixi"

    ```bash
    pixi install -c bioconda kaptive
    ```


## :lucide-blocks: Dependencies

As of `v3.3.0`, Kaptive does not rely on any external binaries and is 100% self-contained.
We do rely on the following dependencies, as defined in the `pyproject.toml`:

 - [python](https://www.python.org/): The core language of the library.
 - [numpy](https://numpy.org/): For vectorised numerical operations.
 - [numba](https://numba.pydata.org/): For speeding up heavy algebra.
 - [rammappy](https://tomdstanton.github.io/rammappy/): For fast, sensitive, nucleotide alignment.
 - [gb-io](https://gb-io.readthedocs.io/): For parsing databases.

### :lucide-package: Optional Dependencies

Kaptive comes with optional dependency groups to keep installation light for
most users. These can be easily installed with `(uv|pip) install kaptive[group1,group2,...]`.


=== "plot"
    
    This installation comes with [plotly](https://plotly.com/python/) for interactive visualisation of results.

    ```bash
    pip install kaptive[plot]
    ```

=== "json"

    This installation comes with [orjson](https://github.com/ijl/orjson) for fast serialsation/deserialsation of results containing numpy arrays. It's mostly used by Kaptive-Web but it's needed in the CLI when dumping serialised results to a JSONL file.

    ```bash
    pip install kaptive[json]
    ```

