---
title: Overview
author: Tom Stanton
comments: true
tags: [databases]
icon: lucide/database
categories:
  - Databases
---

<a id="locus-definition"></a>

## :lucide-help-circle: What is a locus?

A locus in the Kaptive sense refers to a biosynthetic gene cluster that
is responsible for the synthesis of a bacterial surface polysaccharide,
e.g. the *Klebsiella pneumoniae* K locus is responsible for the
synthesis of the capsular polysaccharide, also known as the K antigen.
Each locus in the Kaptive databases has been defined based on a unique
set of genes, with the assumption that this encodes a unique
polysaccharide structure. In many cases, these unique structures will
result in unique immunological serotypes.

The gene translations (protein sequences) from each locus are compared
by pairwise alignment, and must fall under a defined percent identity
threshold to be considered 'unique'. Some genes (such as the core
assembly machinery) will be highly similar, however the genes
responsible for the polysaccharide structural diversity are expected to
be more variable. The default identity thresholds vary across species and are defined by the database curators within the database [metadata file](./curation.md).


## :lucide-database: Available Databases

Below is a list of Kaptive databases you can install using their keyword:

| Species | Locus | Keyword |
|:-:|:-:|:-:|
| [_Klebsiella pneumoniae_ Species Complex](https://github.com/klebgenomics/KpSC_surface_antigen_loci) | K | `kpsc_k` |
| [_Klebsiella pneumoniae_ Species Complex](https://github.com/klebgenomics/KpSC_surface_antigen_loci) | O | `kpsc_o` |
| [_Klebsiella oxytoca_ Species Complex](https://github.com/klebgenomics/KoSC-surface-antigen-loci) | K | `kosc_k` |
| [_Klebsiella oxytoca_ Species Complex](https://github.com/klebgenomics/KoSC-surface-antigen-loci) | O | `kosc_o` |
| [_Acinetobacter baumannii_](https://github.com/johannajkenyon/Abaumannii_surface_polysaccharide_loci) | K | `ab_k` |
| [_Acinetobacter baumannii_](https://github.com/johannajkenyon/Abaumannii_surface_polysaccharide_loci) | OC | `ab_o` |
| [_Escherichia coli_](https://github.com/rgladstone/EC-K-typing) | Group 2 + 3 CPS | `ecoli_kps` |

### :lucide-terminal: Managing Databases via the CLI

Kaptive provides a dedicated command group, [`kaptive db`](../cli/db.md), for managing your local databases.

??? example "Install a Known Database"
    === "Specific database"
        To install an officially supported database, use its known keyword:
        ```bash
        kaptive db install kpsc_k
        ```
    === "All databases"
        You can also install all officially supported databases concurrently:
        ```bash
        kaptive db install all
        ```

??? example "Add a Custom Database"

    To add a new custom database hosted on GitHub, you need the repository owner, repository name, and the base name of the database files.

    For example, if your database files are `My_Custom_Loci.gbk` and `My_Custom_Loci.toml` in `my-org/my-kaptive-db`:
    ```bash
    kaptive db add My_Custom_Loci my-org my-kaptive-db
    ```

??? example "Update Databases"
    To check for and install updates for all installed databases:
    === "Specific database"
        ```bash
        kaptive db update kpsc_k
        ```
    === "All databases"
        ```bash
        kaptive db update
        ```

??? example "List Installed Databases"
    To see which databases are currently installed in your `~/.kaptive` cache:
    ```bash
    kaptive db ls
    ```

??? example "Reset Cache"
    If you want to free up space or start fresh, you can uninstall all local databases and clear the cache:
    ```bash
    kaptive db reset
    ```

---

### :lucide-code: Managing Databases via the Python API

If you are using Kaptive programmatically as a Python library, you can interact with the database system using the 
[`Database`][kaptive.db.Database] or [`DatabaseManager`][kaptive.db.DatabaseManager] classes from
the [`kaptive.db`][kaptive.db] module.

??? example "Installing and Loading"
    ```python
    from kaptive.db import Database, DatabaseManager

    db = Database.install("kpsc_k")  # (1)!
    dbs = Database.install("all")  # (2)!
    db = Database.load("kpsc_k")  # (3)!
    ```

    1. Install a known database (compiles and caches it)
    2. Install all known databases concurrently
    3. Load an already installed database from the local cache

??? example "Adding a Custom Database"
    You can fetch and compile a custom database directly from a GitHub repository:
    ```python
    db = Database.add(
        owner="klebgenomics",
        repo_name="KpSC_surface_antigen_loci",
        db_name="Klebsiella_pneumoniae_Species_Complex_K",
        branch="main",
    )
    ```

??? example "Updating"
    ```python
    updated_dbs = DatabaseManager.update("all")  # (1)!
    updated_db = DatabaseManager.update("kpsc_k")  # (2)!
    ```

    1. Check and update all installed databases
    2. Update a specific database
