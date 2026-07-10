# Kaptive Databases

<a id="locus-definition"></a>

## What is a locus?

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
be more variable. The specific identity thresholds vary across species.

## Format

### Genbank file

Kaptive stores databases in Genbank format consisting of **unique** loci
each with a single record with the following requirements:

- The `source` feature must contain a `note` qualifier which begins with
  a label such as `K locus:`. Whatever follows is used as the locus name
  reported in the Kaptive output. The label is automatically determined,
  and any consistent label ending in a colon will work. However, the
  user can specify exactly which label to use with `--locus_label`, if
  desired.
- The `source` feature may optionally contain a `note` qualifier which
  begins with a label such as `K type:` that specifies the
  phenotype (e.g. serotype or defined polysaccharide structure name) associated with the locus (if known). In cases where only
  some loci are associated with known phenotypes we recommend adding a
  `note` such as `K type: unknown`. If no `type` notes are specified for
  any loci, Kaptive will list them as `unknown` in the output.
  (Kaptive v2.0+)
- Any locus gene should be annotated as `CDS` features. All `CDS`
  features will be used and any other type of feature will be ignored.
- If the gene has a name, it should be specified in a `gene` qualifier.
  This is not required for Kaptive to run, but if absent the gene will
  only be named using its numbered position in the locus and it will not
  be checked for any specific sequence variations relevant to [phenotype
  prediction](Databases.md#phenotype-logic).

Example piece of input Genbank file:

    source          1..23877
                    /organism="Klebsiella pneumoniae"
                    /mol_type="genomic DNA"
                    /note="K locus: KL1"
                    /note="K type: K1"
    CDS             1..897
                    /gene="galF"

#### Nomenclature

In constructing the databases included with Kaptive, we have used the
following nomenclature rules:

- Loci are named after their respective antigen (**K**, **O**, or
  **OC**) followed by the letter **L** (which stands for **Locus**),
  which separates the label for the genotype from the phenotype (e.g.
  KL1 -\> K1). These letters should be in upper case.
- Loci are numbered, first, by their corresponding phenotype, and second,
  in the order in which they were discovered. For example, *Klebsiella*
  K-loci KL1-KL72, KL74 and KL79-KL82 correspond to the originally defined K-types K1-K72, K74 and K79-K82, respectively. K-loci 101 and greater
  correspond to K-loci for which the phenotypes were unknown at the time of locus discovery, numbered in the order in which they
  were discovered.
- Locus genes are named in three parts delimited by an **underscore**
  (**\_**):
  1.  The locus the gene belongs to, e.g. `KL1_` for a gene in the `KL1`
      locus.
  2.  The position of the gene in the locus, e.g. `KL1_01` for the first
      gene in the `KL1` locus.
  3.  The name of the gene as a three-letter italicized symbol written
      in lower case letters and usually suffixed with an italicized
      capital letter, e.g. `KL1_01_galF` for the *galF* gene in the
      `KL1` locus. If the gene name is unknown, this part will be blank
      and the gene instead would be called `KL1_01`.

!!! note
    Databases **must** follow this nomenclature system for distribution
    within Kaptive.

## Metadata 📀
All Kaptive databases now must be accompanied by a metadata [TOML file](https://toml.io/) with the **same name** as the corresponding Genbank file, 
with a '.toml' extension in place of the '.gbk' extension.

Below is an example of the _Klebsiella pneumoniae_ Species Complex K-locus database metadata:

```toml
name = "Klebsiella_pneumoniae_Species_Complex_K"
keyword = "kpsc_k"
genbank = "Klebsiella_pneumoniae_Species_Complex_K.gbk"
organism = "Klebsiella pneumoniae Species Complex"
taxon = 3390273
antigen = "Capsular polysaccharide"
pathway = "Wzx/Wzy-dependent"
prefix = "K"
version = "3.2.1"
id_threshold = 82.5
doi = ["TBD"]
owner = "klebgenomics"
repo = "KpSC_surface_antigen_loci"
branch = "main"
contact = { "Kelly Wyres" = "kaptive.typing@gmail.com" }

[phenotype_logic]
"Capsule null" = { loci = ["KL*"], inactive_genes = ['wza','wzb','wzc','wzx','wzy', 'wcaJ*', 'wbaP*'], priority = 100 }
"K37" = { loci = ["KL22"], inactive_genes = ["atr12"] }
```

<a id="phenotype-logic"></a>
### Phenotype Logic 🧠
The [TOML format](https://toml.io/) is a simple, human-readable, easily-parsable format, which makes it perfect for metadata. For these reasons, it also made sense to define the
phenotype logic here too! Whilst this is still a work-in-progress, here is how we're currently defining it:

1. Each line represents a unique phenotype that can be applied to a serotyping call.
1. All fields accept a wildcard ([`*`](https://docs.python.org/3/library/fnmatch.html)) for selecting multiple items.
1. Loci are defined by the "loci" field - here you can choose the specific loci the logic applies to.

Phenotype logic (previously called "special logic") is a set of rules
that Kaptive uses to predict the polysaccharide phenotype based on the
genes it finds. This was initially implemented for the *Klebsiella
pneumoniae* O locus, whereby additional genes outside of the locus are
used to predict the O antigen (sub)type. This logic was extended to the
*A. baumannii* K locus in Kaptive v2.0.2.

In Kaptive 3, we thought about how we could extend this given what we
know about truncations or other sequence variations of specific genes in
the locus and the impact on the phenotype. For example, in the
*Klebsiella pneumoniae* K locus, we know that a truncation of the core
initiating glycosyltransferase (*wcaJ*) results in a capsule-null
phenotype.

The relevant sequence variations are detailed in the database logic
files, each labelled with the same file prefix as its respective locus
database, and marked with the extension `.logic`. Each line consists of
three tab-separated columns and represents a phenotype rule:

1.  **loci** - the loci the rule applies to (or *ALL* if the rule
    applies to all loci in the database)
2.  **genes** - the genes (and optional state) the rule applies to (or
    *ALL* if the rule applies to all genes in the locus)
3.  **phenotype** - the resulting phenotype that appears in the
    <span class="title-ref">Type</span> column of the Kaptive tabular
    output, replacing the default phenotype i.e. the one specified in
    the locus genbank source identifier in the matching locus database.

Let's look at an example of a logic file for the *K. pneumoniae* K
locus:

| loci | genes             | phenotype    |
|------|-------------------|--------------|
| ALL  | wcaJ1,truncated   | Capsule null |
| KL22 | KL22_17,truncated | K37          |

In the first line, you can see that if *wcaJ1* is truncated in any locus
(selected with *ALL*), the phenotype will be predicted as 'Capsule
null'. Here, any gene with the name *wcaJ1* will be considered, and the
state of the gene is specified as *truncated*. In the last line, you can
see that if *KL22_17* (acetyl-transferase) is truncated in locus KL22,
the phenotype is predicted as 'K37', the non-acetylated version of the
K22 capsule.

!!! note
    The gene name and state are delimited by a comma.


!!! note
    The default phenotype is the "type" label in the Genbank record (e.g.
    K1).


Let's look at an example that uses extra genes outside of the locus
(from the *K. pneumoniae* O locus database):

| loci                 | genes | phenotype |
|----------------------|-------|-----------|
| OL2α.1;OL2α.2;OL2α.3 | orf8  | O2αγ      |

Here, the first line states that if *orf8* is present in a genome
carrying any of the OL2α.1, OL2α.2 or OL2α.3 loci, the phenotype will be
predicted as 'O2αγ'.

!!! note
    Each specific locus and gene is delimited by a semicolon.


!!! note
    Default state is 'presence'.


This logic is applied during the [phenotype
prediction](Method.md#phenotype-prediction) step of typing and is
reported in the <span class="title-ref">Type</span> column of the
Kaptive tabular output.

### App 💫
We have created a simple [![Streamlit App](https://img.shields.io/badge/Streamlit-%23FE4B4B.svg?logo=streamlit&logoColor=white)](https://kaptive-database-validator.streamlit.app/) App to help you 
generate the metadata any database in your repo!


## Distribution

Kaptive 3 introduces a decentralized database system. Instead of bundling massive database files directly into the Kaptive application source code, databases are now decoupled and hosted in their own independent GitHub repositories. This provides several major advantages:

- **Independent Release Cycles**: Database curators can publish updates to a database without waiting for a new Kaptive application release.
- **Custom Databases**: Users can easily create and maintain their own private or public databases in GitHub and add them directly to their local Kaptive installation.
- **Smaller Application Footprint**: You only download the databases you actually need, keeping Kaptive's base installation lightweight.

### How It Works Behind the Scenes

When you install or update a database, Kaptive performs the following sequence of operations under the hood (managed in `kaptive.db.DatabaseManager`):

1. **Metadata Fetch**: Kaptive connects to the specified GitHub repository (via `raw.githubusercontent.com`) and fetches the small `.toml` metadata file.
2. **Version Comparison**: It compares the `version` defined in the remote TOML file against your local installation cache. If your local version is equal to or greater than the remote version, Kaptive skips the download, saving bandwidth and time.
3. **Database Download**: If a new version is detected, Kaptive downloads the full GenBank (`.gbk`) file.
4. **Compilation**: Kaptive parses the GenBank sequences and compiles them into a highly optimized, flat Structure-of-Arrays (SoA) layout (`Database` object). This vectorised layout drastically reduces cache misses during ML and alignment steps.
5. **Caching**: The compiled database is serialized (pickled) and saved to your local cache directory (`~/.kaptive/<keyword>.pkl`), alongside a lightweight JSON file (`~/.kaptive/<keyword>.json`) used for ultra-fast version checking on subsequent runs.

---

### Managing Databases via the CLI

Kaptive provides a dedicated command group, `kaptive db`, for managing your local databases.

#### Install a Known Database
To install an officially supported database, use its known keyword:
```bash
kaptive db install kpsc_k
```
You can also install all officially supported databases concurrently:
```bash
kaptive db install all
```

#### Add a Custom Database
To add a new custom database hosted on GitHub, you need the repository owner, repository name, and the base name of the database files.
```bash
kaptive db add <db_name> <owner> <repo_name> --branch main
```
For example, if your database files are `My_Custom_Loci.gbk` and `My_Custom_Loci.toml` in `my-org/my-kaptive-db`:
```bash
kaptive db add My_Custom_Loci my-org my-kaptive-db
```

#### Update Databases
To check for and install updates for all installed databases:
```bash
kaptive db update
```
Or for a specific database:
```bash
kaptive db update kpsc_k
```

#### List Installed Databases
To see which databases are currently installed in your `~/.kaptive` cache:
```bash
kaptive db ls
```

#### Reset Cache
If you want to free up space or start fresh, you can uninstall all local databases and clear the cache:
```bash
kaptive db reset
```

---

### Managing Databases via the Python API

If you are using Kaptive programmatically as a Python library, you can interact with the database system using the `Database` or `DatabaseManager` classes from the `kaptive.db` module.

#### Installing and Loading
```python
from kaptive.db import Database, DatabaseManager

# Install a known database (compiles and caches it)
db = Database.install("kpsc_k")

# Install all known databases concurrently
dbs = Database.install("all")

# Load an already installed database from the local cache
db = Database.load("kpsc_k")
```

#### Adding a Custom Database
You can fetch and compile a custom database directly from a GitHub repository:
```python
db = Database.add(
    owner="klebgenomics",
    repo_name="KpSC_surface_antigen_loci",
    db_name="Klebsiella_pneumoniae_Species_Complex_K",
    branch="main"
)
```

#### Updating
```python
# Check and update all installed databases
updated_dbs = DatabaseManager.update("all")

# Update a specific database
updated_db = DatabaseManager.update("kpsc_k")
```


## Database Versioning & Release Workflow 🚀
This repository uses a fully automated Continuous Integration / Continuous Deployment (CI/CD) pipeline to manage database versions.

You do not need to manually edit version numbers or create Git tags. The pipeline relies on Semantic Versioning (SemVer) and reads your 
commit messages to automatically calculate the correct version bump, update the corresponding .toml files, and generate 
database-specific release tags.

### How It Works: Conventional Commits ⚙️
The automation script decides how to version a database based on the language used in your commit messages. 
We follow the [Conventional Commits standard](https://www.conventionalcommits.org/en/v1.0.0/).

When you commit changes to a database Genbank file, prefix your commit message with one of the following:

#### Patch Bump 🔨
`fix:` - Use this for correcting typos, fixing broken logic rules, or minor backwards-compatible bug fixes.

- Example: `fix: correct wcaJ truncation rule in Klebsiella`
- Result: `v3.2.1 ➡️ v3.2.2`

#### Minor Bump 🛠️
`feat:` - Use this when adding new features, such as adding a new locus, a new glycosidic linkage, or expanding the phenotype logic in a backwards-compatible way.

- Example: `feat: add KL102 locus to Klebsiella_pneumoniae_K`
- Result: `v3.2.1 ➡️ v3.3.0`

#### Major Bump 🧰
`feat!:` or `[major]` - Use this for breaking changes, such as overhauling the TOML schema, changing existing core 
nomenclature, or deleting previously supported loci.

- Example: `feat!: restructure TOML schema for phenotype logic`
- Result: `v3.2.1 ➡️ v4.0.0`

#### No Bump 🤷
`chore:`, `docs:`, `style:` - Changes to `README`s, generic repository maintenance, or formatting will not trigger a version bump.

### Day-to-Day Workflows 🖇️
#### Updating an Existing Database ⬆️
To update an existing database, simply make your changes to the .gbk files and commit them using the appropriate prefix.

```bash
# 1. Make changes to your files
git add Klebsiella_pneumoniae_K.gbk Klebsiella_pneumoniae_K.toml

# 2. Commit using a Conventional Commit message
git commit -m "feat: add new Wzy-dependent linkage rules"

# 3. Push to main
git push origin main
```

**What happens next?** The GitHub Action will detect the changes to the Klebsiella files, parse the `feat:` prefix, 
bump the minor version in `Klebsiella_pneumoniae_K.toml`, commit that TOML update back to the repository, 
and create a scoped tag (e.g., `Klebsiella_pneumoniae_K-v3.3.0`).

#### Adding a Completely New Database ➡️
The pipeline is database-agnostic. To add a new database, you just need to drop the required files into the repository.

1. Add your new GenBank file (e.g., `seudomonas_aeruginosa_O.gbk`).
2. Add a starting TOML file (e.g., `Pseudomonas_aeruginosa_O.toml`) and manually set the initial version (e.g., `version = "1.0.0"`).
3. Commit and push:

```bash
git add Pseudomonas_aeruginosa_O.*
git commit -m "feat: initial release of Pseudomonas O-locus database"
git push origin main
```

The pipeline will automatically discover the new `.toml` file, register the `feat:` bump (e.g., `v1.1.0`), and tag it.

#### Updating Multiple Databases at Once ⬆️⬆️
If you make a broad change that affects multiple databases (for example, fixing a shared logic rule across both Klebsiella_pneumoniae_K and Klebsiella_pneumoniae_O), simply commit them together:

```bash
git add *.logic
git commit -m "fix: standardize capsule null logic across all databases"
git push origin main**
```

The workflow will detect every database that was modified, bump their `.toml` versions independently, and generate a separate release tag for each one.

### Important Rules ⚠️
 - Never manually edit the version = "..." string in the .toml files. The Python automation (tomlkit) handles this to ensure
   strict alignment between the file contents and the Git tags.
 - Ensure file names match exactly. The base name of the TOML file must match the base name of the GenBank files
   (e.g., `Database_Name.toml` pairs with `Database_Name.gbk`).
 - Pull before you work. Because the GitHub Action makes automated commits to update the TOML files, always
   run `git pull` before starting new work to ensure your local branch has the latest version strings.

## References 📚
[^1]: Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and
    Accurate in silico Antigen Typing with Kaptive 3.
    2025 _Microbial Genomics_ 11(6):001428.
    <https://doi.org/10.1099/mgen.0.001428>
