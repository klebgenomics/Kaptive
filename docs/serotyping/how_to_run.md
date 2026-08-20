---
title: How to run
author: Tom Stanton
comments: true
tags: [serotyping, CLI]
icon: lucide/files
categories:
  - Serotyping
---


## :lucide-file-text: Inputs

Kaptive performs _in silico_ serotyping on bacterial whole genome assemblies. Your input genome assemblies should be in FASTA format (they can be gzip-compressed), with one assembly per file.

```
  >contig1
  ATGAAAAAATGCGCGTATGACT
  >contig2
  TTCGACTCCTGACTGACTGACTTTTATTT
```

## :lucide-dna: Go serotyping!

Once you have [installed Kaptive](../index.md#1-install-kaptive) and your chosen database, you are ready to serotype.

To run Kaptive on your assemblies you can run: 

```bash
kaptive type kpsc_k -o results.tsv *.fasta 
```

This assumes your input files have the `.fasta` file extenion and are avialable in the current directory.  
Kaptive will run on each assembly and print the output to a tab-delimted file called `results.tsv`.    
`kpsc_k` is the database keyword, and points to the _Klebsiella pneumoniae_ Species Complex K database. For a full list of supported key words see [here](../db/overview.md#available-databases). You can also run Kaptive on a custom database by following [these instructions](../db/overview.md#managing-databases-via-the-cli). 

## :lucide-files: Outputs

### :lucide-table: Tabular

The main output of the assembly typing mode is a tab-delimited table of
the results. See [here](./results.md) for tips on interpreting these results. For full explanation of the column content see [here][kaptive.serotyping.io].

The default is to print this table to **stdout**. You can use UNIX
redirection operators (`>` or `>>`) or the `-o`/`--out` flag to write to
a file.

If the summary table already exists and is not empty, Kaptive will
append to it (not overwrite it) and suppress the header line. This
allows you to run Kaptive in succession on sets of assemblies, all
outputting to the same table file.

To disable the tabular output, simply redirect the output to
`/dev/null`.

<a id="fasta"></a>

### :lucide-file-text: Locus sequences

The `-l/--loci` flag produces a fasta file of the region(s) of the
assembly which correspond to the best locus match. This may be a single
piece (in cases of a good assembly and a strong match) or it may be in
multiple pieces (in cases of poor assembly and/or a novel locus).

You can specify either a directory, which will write one file per
assembly named as `{assembly}_kaptive_results.fna`, **or** a single file
("-" for `stdout`), which will write all the sequences to that file.

For example:

```bash
    kaptive type kpsc_k assembly.fasta -l
```

This results in default behaviour which will produce one file per
assembly in the current directory. However, to specify a directory:

```bash
    kaptive type kpsc_k assembly.fasta -l kaptive_results/
```
or for a single file, both are valid:

```bash
    kaptive type kpsc_k assembly.fasta -l kaptive_results.fna
    kaptive type kpsc_k assembly.fasta -l - > kaptive_results.fna
```

!!! note
    This is the same as the `--fna` flag in `kaptive convert`.


<a id="json"></a>

### :lucide-arrow-big-right: Locus plots

Kaptive can generate interative plots showing the locus pieces and genes present in your input assemblies.

```bash
    kaptive type kpsc_k *.fasta --plots my_plot_directory
```

If no directory is specified, Kaptive will generate the plots in the current directory.

!!! note
    Make sure you have installed the correct [dependencies](../index.md#optional-dependencies) to support plotting. 


### :lucide-notebook-tabs: PHA4GE genotyping spec

The [Public Health Alliance for Genomic Epidemiology](https://pha4ge.org/) has developed the [PHA4GE Microbial Genotyping Data Specification](https://github.com/pha4ge/genotyping-specification#about), which represents a standardised format for communicating genotyping methods and results. Kaptive can optionally output results in this format:

```bash
    kaptive type kpsc_k *.fasta --pha4ge my_pha4ge_table.tsv
```

### :lucide-braces: JSON

The `-j/--json` flag produces a JSON file of the results which allows
Kaptive to reconstruct the `TypingResult` objects after a run which can
be used with [kaptive-convert](../cli/convert.md). Unlike
previous version (2 and below), this is a JSON lines file (or "-" for
`stdout`), where each line is a JSON object representing the results for
a single assembly. If the file already exists, Kaptive will append to it
(not overwrite it).

The default is to write this file to: `kaptive_results.json`, however
the path can be specified after the flag, for example:

    kaptive assembly kpsc_k assembly.fasta -j kaptive_results.json

!!! warning
    It is possible to write **all** text formats (TSV, JSON and FASTA) to
    the same file (including stdout), however this is not recommended for
    downstream analysis.
