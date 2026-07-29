---
title: Outputs
author: Tom Stanton
comments: true
tags: [serotyping, CLI]
icon: lucide/files
categories:
  - Serotyping
---

## Tabular

The main output of the assembly typing mode is a tab-delimited table of
the results. See [the api](api.md#kaptive.serotyping.KaptiveRow) for an
explanation of the column content.

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

## Fasta

The `-f/--fasta` flag produces a fasta file of the region(s) of the
assembly which correspond to the best locus match. This may be a single
piece (in cases of a good assembly and a strong match) or it may be in
multiple pieces (in cases of poor assembly and/or a novel locus).

You can specify either a directory, which will write one file per
assembly named as `{assembly}_kaptive_results.fna`, **or** a single file
("-" for `stdout`), which will write all the sequences to that file.

For example:

    kaptive assembly kpsc_k assembly.fasta -f

This results in default behaviour which will produce one file per
assembly in the current directory. However, to specify a directory:

    kaptive assembly kpsc_k assembly.fasta -f kaptive_results/

or for a single file, both are valid:

    kaptive assembly kpsc_k assembly.fasta -f kaptive_results.fna
    kaptive assembly kpsc_k assembly.fasta -f - > kaptive_results.fna

!!! note
    This is the same as the `--fna` flag in `kaptive convert`.


<a id="json"></a>

## JSON

The `-j/--json` flag produces a JSON file of the results which allows
Kaptive to reconstruct the `TypingResult` objects after a run which can
be used with [kaptive-convert](cli.md#kaptive-convert). Unlike
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
