---
title: database (db)
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/database
categories:
  - Databases
---

📦 Manage local and remote Kaptive databases.

Aggregates subcommands for listing, installing, updating, resetting, adding,
extracting, and displaying metadata for Kaptive databases.

Aliases:
    db


```text
[1;36musage: [0m[1;36m[0mkaptive database [1;36m[options][0m [1;36m[subcommand][0m ...

[1m📦 Manage local and remote Kaptive databases.

Aggregates subcommands for listing, installing, updating, resetting, adding,
extracting, and displaying metadata for Kaptive databases.

Aliases:
    db
[0m

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress

[1;36m[1m'database' subcommands[0m[0m:
    list (ls)           📋 List all currently installed local databases.
    available (avail)   🌐 List all available official databases for installation.
    add                 🔗 Add a custom reference database from a GitHub repository.
    install             📦 Install known reference databases via keyword.
    update              🔄 Update installed local databases from remote repositories.
    reset               🧹 Uninstall all local databases and reset local cache.
    extract             📤 Extract database records in FASTA format.
    metadata (info)     📊 Print detailed metadata of a Kaptive database.
```

## Subcommands

### list (ls)

📋 List all currently installed local databases.

Displays the keywords of all compiled `.pkl` databases found in the user's
local Kaptive directory (`~/.kaptive`).

Aliases:
    ls


```text
[1;36musage: [0m[1;36m[0m[1;36m[0mkaptive database list [1;36m[options][0m

[1m📋 List all currently installed local databases.

Displays the keywords of all compiled `.pkl` databases found in the user's
local Kaptive directory (`~/.kaptive`).

Aliases:
    ls
[0m

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### available (avail)

🌐 List all available official databases for installation.

Displays the keywords of all officially supported databases curated in GitHub repositories
that can be installed via `kaptive db install <keyword>`.

Aliases:
    avail


```text
[1;36musage: [0m[1;36m[0m[1;36m[0mkaptive database available [1;36m[options][0m

[1m🌐 List all available official databases for installation.

Displays the keywords of all officially supported databases curated in GitHub repositories
that can be installed via `kaptive db install <keyword>`.

Aliases:
    avail
[0m

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### add

🔗 Add a custom reference database from a GitHub repository.

Fetches GenBank and TOML metadata files from any specified GitHub owner/repo/branch,
compiles the database, and registers it in the local cache.


```text
[1;36musage: [0m[1;36m[0m[1;36m[0mkaptive database add [1;36m[options][0m database owner repo_name

[1m🔗 Add a custom reference database from a GitHub repository.

Fetches GenBank and TOML metadata files from any specified GitHub owner/repo/branch,
compiles the database, and registers it in the local cache.
[0m

[1;36m📥 Inputs[0m:
  database              Name for the new database

[1;36m[1m🌐 GitHub Details[0m[0m:
  owner                 GitHub repository owner
  repo_name             GitHub repository name
  -b, --branch [BRANCH]
                        GitHub repository branch (default: main)

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress
```

### install

📦 Install known reference databases via keyword.

Downloads GenBank and TOML definition files from official repositories, compiles
them into vectorized [`Database`][kaptive.db.core.Database] objects, and caches
them locally.


```text
[1;36musage: [0m[1;36m[0m[1;36m[0mkaptive database install [1;36m[options][0m database

[1m📦 Install known reference databases via keyword.

Downloads GenBank and TOML definition files from official repositories, compiles
them into vectorized [`Database`][kaptive.db.core.Database] objects, and caches
them locally.
[0m

[1;36m📥 Inputs[0m:
  database       Database keyword (see: `kaptive db avail`) or 'all'

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### update

🔄 Update installed local databases from remote repositories.

Checks installed databases against their source GitHub repositories for newer
versions defined in TOML metadata and re-compiles modified databases.


```text
[1;36musage: [0m[1;36m[0m[1;36m[0mkaptive database update [1;36m[options][0m [database]

[1m🔄 Update installed local databases from remote repositories.

Checks installed databases against their source GitHub repositories for newer
versions defined in TOML metadata and re-compiles modified databases.
[0m

[1;36m📥 Inputs[0m:
  database       Database keyword (see: `kaptive db list`) or 'all' (default: all)

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### reset

🧹 Uninstall all local databases and reset local cache.

Deletes all compiled `.pkl` and metadata `.json` files from the local Kaptive cache directory.


```text
[1;36musage: [0m[1;36m[0m[1;36m[0mkaptive database reset [1;36m[options][0m

[1m🧹 Uninstall all local databases and reset local cache.

Deletes all compiled `.pkl` and metadata `.json` files from the local Kaptive cache directory.
[0m

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### extract

📤 Extract database records in FASTA format.

Aggregates subcommands for extracting locus nucleotide sequences (`loci`), gene
coding sequences (`genes`), and translated amino acid sequences (`proteins`).


```text
[1;36musage: [0m[1;36m[0m[1;36m[0mkaptive database extract [1;36m[options][0m [1;36m[subcommand][0m ...

[1m📤 Extract database records in FASTA format.

Aggregates subcommands for extracting locus nucleotide sequences (`loci`), gene
coding sequences (`genes`), and translated amino acid sequences (`proteins`).
[0m

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress

[1;36m[1m'extract' subcommands[0m[0m:
    loci                🧬 Extract locus nucleotide sequences in FASTA format.
    genes               🧩 Extract gene coding sequences in FASTA format.
    proteins            🧶 Extract translated protein sequences in FASTA format.
```

### metadata (info)

📊 Print detailed metadata of a Kaptive database.

Displays summary information including organism, taxon ID, antigen type, synthesis pathway,
version, identity threshold, GenBank filename, DOIs, repository URL, and curator contacts.

Aliases:
    info


```text
[1;36musage: [0m[1;36m[0m[1;36m[0mkaptive database metadata [1;36m[options][0m database

[1m📊 Print detailed metadata of a Kaptive database.

Displays summary information including organism, taxon ID, antigen type, synthesis pathway,
version, identity threshold, GenBank filename, DOIs, repository URL, and curator contacts.

Aliases:
    info
[0m

[1;36m📥 Inputs[0m:
  database       Database path or keyword (see: `kaptive db list`)

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```
