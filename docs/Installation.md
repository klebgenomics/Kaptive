# Installation

## Dependencies

As of v3.3, Kaptive does not rely on any external binaries and is 100% self-contained.
We do rely on the following dependencies, as defined in the `pyproject.toml`:
 - [python](https://www.python.org/): The core language of the library.
 - [numpy](https://numpy.org/): For vectorised numerical operations.
 - [numba](https://numba.pydata.org/): For speeding up heavy algebra.
 - [mappy](https://github.com/lh3/minimap2/tree/master/python): For mapping genes to contigs.
 - [gb-io](https://gb-io.readthedocs.io/): For parsing external databases.

## Download and install Kaptive

### From [PyPI](https://pypi.org/project/kaptive/):

#### via [pip](https://pypi.org/project/pip/):
```bash
$ pip install kaptive
```

#### via [uv](https://docs.astral.sh/uv/):
```bash
$ uv pip install kaptive
# OR
$ uv add kaptive
```

### From [Bioconda](https://anaconda.org/bioconda/kaptive)

#### via [conda]() or [mamba]():
```bash
$ conda install bioconda::kaptive
# OR
$ mamba install bioconda::kaptive
```

#### via `pixi`:
```bash
$ pixi install -c bioconda kaptive
```

## Optional Dependencies

Kaptive comes with optional dependency groups to keep installation light for
most users. These can be easily installed with `(uv|pip) install kaptive[group1,group2,...]`.

### `kaptive[plot]`
This installation comes with [plotly](https://plotly.com/python/) for interactive
visualisation of results.

### `kaptive[json]`

This installation comes with [orjson](https://github.com/ijl/orjson) for fast serialsation/deserialsation of results containing numpy arrays.
It's mostly used by Kaptive-Web but it's needed in the CLI when dumping serialised results to a JSONL file.

### `kaptive[web]`

This is simply `kaptive[plot,json]`.