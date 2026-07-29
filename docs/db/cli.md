---
title: database (db)
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/terminal
categories:
  - Databases
---

📦 Manage Kaptive databases

```text
usage: kaptive database [options]
                        [subcommand] ...

📦 Manage Kaptive databases

🌎 Global options:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress

'database' subcommands:
    list (ls)           📋 Lists installed databases
    add                 🔗 Add a new database from a GitHub repository.
    install             📦 Install a known database via keyword (or 'all').
    update              🔄 Update installed databases from their remote repositories.
    reset               🧹 Uninstall all local databases and clear the cache.
    extract             📤 Extract records from a known database.
    metadata (info)     📊 Print metadata of a database
```

## Subcommands

### list (ls)

📋 Lists installed databases

```text
usage: kaptive database list [options]

📋 Lists installed databases

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### add

🔗 Add a new database from a GitHub repository.

```text
usage: kaptive database add [options] database owner repo_name

🔗 Add a new database from a GitHub repository.

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

📦 Install a known database via keyword (or 'all').

```text
usage: kaptive database install [options] database

📦 Install a known database via keyword (or 'all').

📥 Inputs:
  database       Database keyword (see: `kaptive db list`) or 'all'

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### update

🔄 Update installed databases from their remote repositories.

```text
usage: kaptive database update [options] [database]

🔄 Update installed databases from their remote repositories.

📥 Inputs:
  database       Database keyword (see: `kaptive db list`) or 'all' (default: all)

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### reset

🧹 Uninstall all local databases and clear the cache.

```text
usage: kaptive database reset [options]

🧹 Uninstall all local databases and clear the cache.

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```

### extract

📤 Extract records from a known database.

```text
usage: kaptive database extract [options] [subcommand] ...

📤 Extract records from a known database.

🌎 Global options:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress

'extract' subcommands:
    loci                🧬 Extracts locus sequences from the database in fasta format.
    genes               🧩 Extracts gene sequences from the database in fasta format.
    proteins            🧶 Extracts protein sequences from the database in fasta format.
```

### metadata (info)

📊 Print metadata of a database

```text
usage: kaptive database metadata [options] database

📊 Print metadata of a database

📥 Inputs:
  database       Database path or keyword (see: `kaptive db list`)

🌎 Global options:
  -h, --help     show this help message and exit
  -V, --verbose  Enable verbose output/progress
```
