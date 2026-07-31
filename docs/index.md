---
title: Home
description: The tool for _in silico_ serotyping
author: Tom Stanton
comments: true
icon: lucide/home
---

<div align="center">
  <img src="docs/assets/logo.png" alt="Kaptive" width="400">
  <p>The tool for <em>in silico</em> serotyping</p>
  <a href="https://pypi.org/project/kaptive/">
    <img src="https://img.shields.io/pypi/pyversions/kaptive" alt="PyPI - Python Version">
  </a>&nbsp;
  <a href="https://anaconda.org/bioconda/kaptive">
    <img src="https://img.shields.io/conda/vn/bioconda/kaptive" alt="Conda Version">
  </a>&nbsp;
  <a href="https://github.com/astral-sh/ruff">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="ruff">
  </a>&nbsp;
  <a href="https://github.com/astral-sh/ty">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ty/main/assets/badge/v0.json" alt="ty">
  </a>&nbsp;
</div>

---

## :lucide-zap: Quickstart

If you're looking to get up and running with *in silico* serotyping of genome assemblies as fast as possible, you're in the right place. 
We'll skip the heavy algorithmic details here and just focus on getting you your results.

### :lucide-download: 1. Install Kaptive

=== ":simple-pypi: pip"

    ```bash
    pip install kaptive
    ```

=== ":simple-pypi: uv"

    ```bash
    uv add kaptive
    ```

=== ":simple-anaconda: conda"

    ```bash
    conda install bioconda::kaptive
    ```

=== ":simple-anaconda: mamba"

    ```bash
    mamba install bioconda::kaptive
    ```

=== ":lucide-box: pixi"

    ```bash
    pixi add kaptive
    ```


### :lucide-database: 2. Download a Database

Kaptive uses beautifully curated databases containing reference surface antigen loci for different species. To see what's available for download, just run:

```bash
kaptive db available
```

To pull down a specific database, use the `install` command with the database name. For example, if you're working with the *Acinetobacter baumannii* K-locus, you'd want `ab_k`:

```bash
kaptive db install ab_k
```

### :lucide-dna: 3. Serotype Your Genomes!

Once Kaptive and your database are installed, you're ready to roll. Your input genome assemblies should be in FASTA format (and don't worry, gzip-compressed files work perfectly too). 

Let's run Kaptive and save the output to a file called `results.tsv`:

```bash
kaptive type ab_k *.fasta > results.tsv
```

### :lucide-bar-chart: 4. Understanding the Output

Kaptive spits out a tab-separated values (TSV) report, which you can easily open up in Excel, Numbers, or any text editor to browse through. 

Here are the most critical columns to keep an eye on in your `results.tsv` file:

* **Assembly**: The name of your input genome file.
* **Best match locus**: The best-matching serotype or locus found in the database (e.g., `KL1`).
* **Confidence**: How confident Kaptive is in this call - this is either "Typeable" or "Untypeable"

For a deeper dive into all the other columns and more advanced outputs, check out the full [Outputs documentation](serotyping/outputs.md). Happy serotyping!


## :lucide-graduation-cap: Tutorial

Step-by-step [video](https://klebnet.org/training/) and
[documented](https://docs.google.com/document/d/1aggwBCGu1CfsduOoKI0e6TRYOGtwwSceSBdKKkjCisA/edit?usp=sharing)
tutorials are available, covering:

- Kaptive's features and their scientific rationale
- How to run Kaptive
- Examples, illustrating how to run and interpret results
- Further investigations (e.g. exploring novel loci, IS insertions)

??? question "Kaptive theory with Kelly!"
    
    <iframe width="560" height="315" src="https://www.youtube.com/embed/aKySfea7ukE?si=6NDD1psd7yDbjOdC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>


??? question "Kaptive usage with Tom!"

    <iframe width="560" height="315" src="https://www.youtube.com/embed/E-ko2QlON8k?si=IyKWSbDnptvhzt6w" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>


!!! note

    The tutorials are based on Kaptive 2.0, but the principles are similar for Kaptive 3.0.


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


=== ":lucide-pencil-sparkles: plot"
    
    This installation comes with [plotly](https://plotly.com/python/) for interactive visualisation of results.

    ```bash
    pip install kaptive[plot]
    ```

=== ":lucide-braces: json"

    This installation comes with [orjson](https://github.com/ijl/orjson) for fast serialsation/deserialsation of results containing numpy arrays. It's mostly used by Kaptive-Web but it's needed in the CLI when dumping serialised results to a JSONL file.

    ```bash
    pip install kaptive[json]
    ```

## :lucide-book-open-text: Citation

!!! abstract "If you use **Kaptive** in your work, please cite:"

    ```bibtex
    @article{mbs:/content/journal/mgen/10.1099/mgen.0.001428,
       author = "Stanton, Thomas David and Hetland, Marit A.K. and Löhr, Iren H. and Holt, Kathryn E. and Wyres, Kelly L.",
       title = "Fast and accurate in silico antigen typing with Kaptive 3", 
       journal= "Microbial Genomics",
       year = "2025",
       volume = "11",
       number = "6",
       pages = "",
       doi = "https://doi.org/10.1099/mgen.0.001428",
       url = "https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001428",
       publisher = "Microbiology Society",
       issn = "2057-5858",
       type = "Journal Article",
       keywords = "serotyping",
       keywords = "Klebsiella",
       keywords = "tools",
       keywords = "capsule",
       keywords = "sero-epidemiology",
       keywords = "antigen",
        eid = "001428",
       abstract = "Surface polysaccharides are common antigens in priority pathogens and therefore attractive targets for novel control strategies such as vaccines, monoclonal antibody and phage therapies. Distinct serotypes correspond to diverse polysaccharide structures that are encoded by distinct biosynthesis gene clusters; e.g. the Klebsiella pneumoniae species complex (KpSC) K- and O-loci encode the synthesis machinery for the capsule (K) and outer-lipopolysaccharides (O), respectively. We previously presented Kaptive and Kaptive 2, programmes to identify K- and O-loci directly from KpSC genome assemblies (later adapted for Acinetobacter baumannii), enabling sero-epidemiological analyses to guide vaccine and phage therapy development. However, for some KpSC genome collections, Kaptive (v≤2) was unable to type a high proportion of K-loci. Here, we identify the cause of this issue as assembly fragmentation and present a new version of Kaptive (v3) to circumvent this problem, reduce processing times and simplify output interpretation. We compared the performance of Kaptive v2 and Kaptive v3 for typing genome assemblies generated from subsampled Illumina read sets (decrements of 10× depth), for which a corresponding high-quality completed genome was also available to determine the ‘true’ loci (n=549 KpSC, n=198 A. baumannii). Both versions of Kaptive showed high rates of agreement to the matched true locus amongst ‘typeable’ locus calls (≥96% for ≥20× read depth), but Kaptive v3 was more sensitive, particularly for low-depth assemblies (at &lt;40× depth, v3 ranged 0.85–1 vs v2 0.09–0.94) and/or typing KpSC K-loci (e.g. 0.97 vs 0.82 for non-subsampled assemblies). Overall, Kaptive v3 was also associated with a higher rate of optimal outcomes; i.e. loci matching those in the reference database were correctly typed, and genuine novel loci were reported as untypeable (73–98% for v3 vs 7–77% for v2 for KpSC K-loci). Kaptive v3 was &gt;1 order of magnitude faster than Kaptive v2, making it easy to analyse thousands of assemblies on a desktop computer, facilitating broadly accessible in silico serotyping that is both accurate and sensitive. The Kaptive v3 source code is freely available on GitHub (https://github.com/klebgenomics/Kaptive), and has been implemented in Kaptive Web (https://kaptive-web.erc.monash.edu/).",
      }
    ```
