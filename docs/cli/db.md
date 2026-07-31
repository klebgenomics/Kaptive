---
title: database (db)
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/database
categories:
  - Databases
---

📦 Manage local and remote Kaptive databases

Aggregates subcommands for listing, installing, updating, resetting, adding,
extracting, and displaying metadata for Kaptive databases.

Aliases:
    db


```text
usage: kaptive database [options] [subcommand] ...

📦 Manage local and remote Kaptive databases

Aggregates subcommands for listing, installing, updating, resetting, adding,
extracting, and displaying metadata for Kaptive databases.

Aliases:
    db

🌎 Global options:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress

'database' subcommands:
    list (ls)           📋 List all currently installed local databases
    available (avail)   🌐 List all available official databases for installation
    add                 🔗 Add a custom reference database from a GitHub repository
    install             📦 Install known reference databases via keyword
    update              🔄 Update installed local databases from remote repositories
    reset               🧹 Uninstall all local databases and reset local cache
    extract             📤 Extract database records in FASTA format
    metadata (info)     📊 Print detailed metadata of a Kaptive database
```

## Subcommands

### list (ls)

📋 List all currently installed local databases

Displays the keywords of all compiled `.pkl` databases found in the user's
local Kaptive directory (`~/.kaptive`).

Aliases:
    ls


```text
usage: kaptive database list [options]

📋 List all currently installed local databases

Displays the keywords of all compiled `.pkl` databases found in the user's
local Kaptive directory (`~/.kaptive`).

Aliases:
    ls

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### available (avail)

🌐 List all available official databases for installation

Displays the keywords of all officially supported databases curated in GitHub repositories
that can be installed via `kaptive db install <keyword>`.

Aliases:
    avail


```text
usage: kaptive database available [options]

🌐 List all available official databases for installation

Displays the keywords of all officially supported databases curated in GitHub repositories
that can be installed via `kaptive db install <keyword>`.

Aliases:
    avail

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### add

🔗 Add a custom reference database from a GitHub repository

Fetches GenBank and TOML metadata files from any specified GitHub owner/repo/branch,
compiles the database, and registers it in the local cache.


```text
usage: kaptive database add [options] database owner repo_name

🔗 Add a custom reference database from a GitHub repository

Fetches GenBank and TOML metadata files from any specified GitHub owner/repo/branch,
compiles the database, and registers it in the local cache.

📥 Inputs:
  database              Name for the new database

🌐 GitHub Details:
  owner                 GitHub repository owner
  repo_name             GitHub repository name
  -b, --branch [BRANCH]
                        GitHub repository branch (default: main)

🌎 Global options:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress
```

### install

📦 Install known reference databases via keyword

Downloads GenBank and TOML definition files from official repositories, compiles
them into vectorized [`Database`][kaptive.db.core.Database] objects, and caches
them locally.


```text
usage: kaptive database install [options] database

📦 Install known reference databases via keyword

Downloads GenBank and TOML definition files from official repositories, compiles
them into vectorized [`Database`][kaptive.db.core.Database] objects, and caches
them locally.

📥 Inputs:
  database       Database keyword (see: `kaptive db avail`) or 'all'

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### update

🔄 Update installed local databases from remote repositories

Checks installed databases against their source GitHub repositories for newer
versions defined in TOML metadata and re-compiles modified databases.


```text
usage: kaptive database update [options] [database]

🔄 Update installed local databases from remote repositories

Checks installed databases against their source GitHub repositories for newer
versions defined in TOML metadata and re-compiles modified databases.

📥 Inputs:
  database       Database keyword (see: `kaptive db list`) or 'all' (default: all)

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### reset

🧹 Uninstall all local databases and reset local cache

Deletes all compiled `.pkl` and metadata `.json` files from the local Kaptive cache directory.


```text
usage: kaptive database reset [options]

🧹 Uninstall all local databases and reset local cache

Deletes all compiled `.pkl` and metadata `.json` files from the local Kaptive cache directory.

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### extract

📤 Extract database records in FASTA format

Aggregates subcommands for extracting locus nucleotide sequences (`loci`), gene
coding sequences (`genes`), and translated amino acid sequences (`proteins`).


```text
usage: kaptive database extract [options] [subcommand] ...

📤 Extract database records in FASTA format

Aggregates subcommands for extracting locus nucleotide sequences (`loci`), gene
coding sequences (`genes`), and translated amino acid sequences (`proteins`).

🌎 Global options:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress

'extract' subcommands:
    loci                🧬 Extract locus nucleotide sequences in FASTA format
    genes               🧩 Extract gene coding sequences in FASTA format
    proteins            🧶 Extract translated protein sequences in FASTA format
```

### metadata (info)

📊 Print detailed metadata of a Kaptive database

Displays summary information including organism, taxon ID, antigen type, synthesis pathway,
version, identity threshold, GenBank filename, DOIs, repository URL, and curator contacts.

Aliases:
    info


```text
usage: kaptive database metadata [options] database

📊 Print detailed metadata of a Kaptive database

Displays summary information including organism, taxon ID, antigen type, synthesis pathway,
version, identity threshold, GenBank filename, DOIs, repository URL, and curator contacts.

Aliases:
    info

📥 Inputs:
  database       Database path or keyword (see: `kaptive db list`)

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```
