---
title: Curation
author: Tom Stanton
comments: true
tags: [databases]
icon: lucide/landmark
categories:
  - Databases
---

## :lucide-tag: Nomenclature

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

!!! danger
    Databases **must** follow this nomenclature system for integration with Kaptive.

## :lucide-book: Reference Loci

Kaptive databases are curated in Genbank format consisting of **unique** loci
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
  prediction](#phenotype-logic).

??? example "Example GenBank snippet"
    
    ```text
    source          1..23877
                    /organism="Klebsiella pneumoniae"
                    /mol_type="genomic DNA"
                    /note="K locus: KL1"
                    /note="K type: K1"
    CDS             1..897
                    /gene="galF"
    ```

## :lucide-file-json: Metadata

All Kaptive databases now must be accompanied by a metadata [TOML file](https://toml.io/) with the **same name** as the corresponding Genbank file, 
with a `.toml` extension in place of the `.gbk` extension.

??? example "Example Metadata TOML"
    An explanation of each metadata field can be found in the [API reference][kaptive.db.DatabaseMetadata].
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
    doi = ["10.1099/mgen.0.000102", "10.1099/mgen.0.000800"]
    owner = "klebgenomics"
    repo = "KpSC_surface_antigen_loci"
    branch = "main"
    contact = { "Wyres lab, wyreslab.com" = "kaptive.typing@gmail.com" }
    ```

### :lucide-network: Phenotype Logic

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

The [TOML format](https://toml.io/) is a simple, human-readable, easily-parsable format, which makes it perfect for metadata.
For these reasons, it also made sense to define the phenotype logic here too!

??? example "Example Inactive Gene Logic"

    ```toml
    [phenotype_logic]
    "Capsule null" = { loci = ["KL*"], inactive_genes = ['wza','wzb','wzc','wzx','wzy', 'wcaJ*', 'wbaP*'], priority = 100 }
    "K37" = { loci = ["KL22"], inactive_genes = ["atr12"] }
    ```

    1. Each line represents a unique phenotype that can be applied to a serotyping call.
    1. All fields accept a wildcard ([`*`](https://docs.python.org/3/library/fnmatch.html)) for selecting multiple items.
    1. Loci are defined by the "loci" field - here you can choose the specific loci the logic applies to.

    In the first line, you can see that if *wcaJ1* is truncated in any locus
    (selected with *ALL*), the phenotype will be predicted as 'Capsule
    null'. Here, any gene with the name *wcaJ1* will be considered, and the
    state of the gene is specified as *truncated*. In the last line, you can
    see that if *KL22_17* (acetyl-transferase) is truncated in locus KL22,
    the phenotype is predicted as 'K37', the non-acetylated version of the
    K22 capsule.

??? example "Example Extra Gene Logic"

    ```toml
    [phenotype_logic]
    "O2β" = { loci = ["OL2α*"], extra_genes = ["gmlA", "gmlB", "gmlC"] }
    "O2αγ" = { loci = ["OL2α*"], extra_genes = ["orf8"] }
    "O1α,2α" = { loci = ["OL2α*"], extra_genes = ["wbbY"] }
    "O1αβ,2α" = { loci = ["OL2α*"], extra_genes = ["wbbY", "wbbZ"] }
    "O1αβ,2γ" = { loci = ["OL2α*"], extra_genes = ["orf8", "wbbY", "wbbZ"] }
    "O1α,2β" = { loci = ["OL2α*"], extra_genes = ["gmlA", "gmlB", "gmlC", "wbbY"] }
    "O1αβ,2β" = { loci = ["OL2α*"], extra_genes = ["gmlA", "gmlB", "gmlC", "wbbY", "wbbZ"] }
    "O11α,2α" = { loci = ["OL2α*"], extra_genes = ["wbmV", "wbmW"] }
    "O11αβ,2α" = { loci = ["OL2α*"], extra_genes = ["wbmV", "wbmW", "wbmX"] }
    "O11αβ,2γ" = { loci = ["OL2α*"], extra_genes = ["orf8", "wbmV", "wbmW", "wbmX"] }
    "O11α,2β" = { loci = ["OL2α*"], extra_genes = ["gmlA", "gmlB", "gmlC", "wbmV", "wbmW"] }
    "O11αβ,2β" = { loci = ["OL2α*"], extra_genes = ["gmlA", "gmlB", "gmlC", "wbmV", "wbmW", "wbmX"] }
    ```

    Here, the first line states that if *gmlA*, *gmlB* and *gmlC* are present in a genome carrying any of the OL2α.1, OL2α.2 or OL2α.3 loci, the phenotype will be predicted as 'O2β'. All genes must be present and detected without truncations in order for this phenotype to apply. If only a subset of genes are detected and/or if one or more of the genes includes a nonsense or frameshift mutation, the phenotype will be reported as the default type specified in the `O type:` label in the O locus `source` qualifier (in this case it would be O2α). 

!!! tip
    We developed an [app](https://kaptive-database-validator.streamlit.app/) to help you validate your database GenBank file, generate the metadata TOML and create phenotype logic.

## :lucide-share-2: Distribution

Kaptive 3 introduces a decentralized database system. Instead of bundling massive database files directly into the Kaptive application source code, databases are now decoupled and hosted in their own independent GitHub repositories. This provides several major advantages:

- **Independent Release Cycles**: Database curators can publish updates to a database without waiting for a new Kaptive application release.
- **Custom Databases**: Users can easily create and maintain their own private or public databases in GitHub and add them directly to their local Kaptive installation.
- **Smaller Application Footprint**: You only download the databases you actually need, keeping Kaptive's base installation lightweight.

### :lucide-settings: How It Works Behind the Scenes

When you install or update a database, Kaptive performs the following sequence of operations under the hood (managed in `kaptive.db.DatabaseManager`):

1. **Metadata Fetch**: Kaptive connects to the specified GitHub repository (via `raw.githubusercontent.com`) and fetches the small `.toml` metadata file.
2. **Version Comparison**: It compares the `version` defined in the remote TOML file against your local installation cache. If your local version is equal to or greater than the remote version, Kaptive skips the download, saving bandwidth and time.
3. **Database Download**: If a new version is detected, Kaptive downloads the full GenBank (`.gbk`) file.
4. **Compilation**: Kaptive parses the GenBank sequences and compiles them into a highly optimized, flat Structure-of-Arrays (SoA) layout (`Database` object). This vectorised layout drastically reduces cache misses during ML and alignment steps.
5. **Caching**: The compiled database is serialized (pickled) and saved to your local cache directory (`~/.kaptive/<keyword>.pkl`), alongside a lightweight JSON file (`~/.kaptive/<keyword>.json`) used for ultra-fast version checking on subsequent runs.


## :lucide-git-branch: Database Versioning & Release Workflow
We have generated a [template Github repository](https://github.com/tomdstanton/kaptive-db-template) that includes a fully automated Continuous Integration / Continuous Deployment (CI/CD) pipeline to manage database versions.  

You do not need to manually edit version numbers or create Git tags. The pipeline relies on Semantic Versioning (SemVer) and reads your 
commit messages to automatically calculate the correct version bump, update the corresponding .toml files, and generate 
database-specific release tags.

The automation script decides how to version a database based on the language used in your commit messages. 
We follow the [Conventional Commits standard](https://www.conventionalcommits.org/en/v1.0.0/):

| Bump | Change | Usage | Prefix | Tag |
|:-:|:-:|:-:|:-:|:-:|
| Major |`x.0.0` | For incompatible or major structural changes. | `major:` | `[major]`
| Minor |`0.x.0` | For adding new features or loci in a backwards-compatible manner. | `minor:` | `[minor]`
| Patch |`0.0.x` | For backwards-compatible bug fixes or minor metadata corrections. | `patch:` | `[patch]`

!!! note
    Any commit missing the above tags will be skipped by the release automation.

### :lucide-square-plus: Establish your database

Setting up your database is easy, just: 

1. Copy the [template repo](https://github.com/tomdstanton/kaptive-db-template)
2. Add your database files  
3. [Customise](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) the `READme` 
(You can see an example customised repo [here](https://github.com/klebgenomics/KpSC_surface_antigen_loci/tree/main).) 

### :lucide-workflow: Example Workflows

??? example "Updating an Existing Database"

    To update an existing database, simply make your changes to the .gbk files and commit them using the appropriate prefix.

    ```bash
    git add Klebsiella_pneumoniae_K.gbk Klebsiella_pneumoniae_K.toml # (1)!
    git commit -m "minor: add new Wzy-dependent linkage rules" # (2)!
    git push origin main # (3)!
    ```

    1. Make changes to your files
    2. Commit using a Conventional Commit message
    3. Push to main

    **What happens next?** The GitHub Action will detect the changes to the Klebsiella files, parse the `minor:` prefix, 
    bump the minor version in `Klebsiella_pneumoniae_K.toml`, commit that TOML update back to the repository, 
    and create a scoped tag (e.g., `Klebsiella_pneumoniae_K-v3.3.0`).

??? example "Adding a Completely New Database"
    
    The pipeline is database-agnostic. To add a new database, you just need to drop the required files into the repository.

    ```bash
    git add Pseudomonas_aeruginosa_O.* # (1)!
    git commit -m "major: initial release of Pseudomonas O-locus database" # (2)!
    git push origin main # (3)!
    ```

    1. Add your new GenBank file (e.g., `Pseudomonas_aeruginosa_O.gbk`).
    2. Add a starting TOML file (e.g., `Pseudomonas_aeruginosa_O.toml`) and manually set the initial version (e.g., `version = "0.0.0"`).
    3. Commit and push!

    The pipeline will automatically discover the new `.toml` file, register the `major:` bump (e.g., `v1.0.0`), and tag it.

??? example "Updating Multiple Databases at Once"
    
    If you make a broad change that affects multiple databases (for example, fixing a shared logic rule across both Klebsiella_pneumoniae_K and Klebsiella_pneumoniae_O), simply commit them together:

    ```bash
    git add *.logic
    git commit -m "fix: standardize capsule null logic across all databases"
    git push origin main**
    ```

    The workflow will detect every database that was modified, bump their `.toml` versions independently, and generate a separate release tag for each one.

!!! warning "Important Rules"
     - Never manually edit the version = "..." string in the .toml files. The Python automation (tomlkit) handles this to ensure
       strict alignment between the file contents and the Git tags.
     - Ensure file names match exactly. The base name of the TOML file must match the base name of the GenBank files
       (e.g., `Database_Name.toml` pairs with `Database_Name.gbk`).
     - Pull before you work. Because the GitHub Action makes automated commits to update the TOML files, always
       run `git pull` before starting new work to ensure your local branch has the latest version strings.
