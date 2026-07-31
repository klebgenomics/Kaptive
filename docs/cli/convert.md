---
title: convert
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/arrow-left-right
categories:
  - Serotyping
---

🔄 Convert serialized Kaptive results into different formats

Reads serialized JSON-lines serotyping output records and converts them into tabular
TSV, PHA4GE TSV, or sequence FASTA files without re-running the serotyping pipeline.


```text
usage: kaptive convert [options] jsonl

🔄 Convert serialized Kaptive results into different formats

Reads serialized JSON-lines serotyping output records and converts them into tabular
TSV, PHA4GE TSV, or sequence FASTA files without re-running the serotyping pipeline.

📥 Inputs:
  jsonl                 Serialised results in JSON-lines format (default: stdin)

📤 Outputs:
  -t, --tsv [FILE]      Write serotyping results as a TSV report to a file (default: stdout)
  -l, --loci [DIR]      Write locus nucleotide fasta files to a directory (default: ./)
  -g, --genes [DIR]     Write gene nucleotide fasta files to a directory (default: ./)
  -p, --proteins [DIR]  Write translation amino-acid fasta files to a directory (default: ./)
  --pha4ge [FILE]       Write PHA4GE-compliant serotyping report to a TSV file (default: kaptive_results.pha4ge)
  --plots [DIR]         Generate interactive locus plots to a directory (default: ./)

🌎 Global options:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress
```
