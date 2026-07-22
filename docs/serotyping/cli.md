---
title: type (assembly)
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/terminal
categories:
  - Serotyping
---

💉 In silico serotyping of assemblies.

```text
usage: kaptive type [options] database genomes [genomes ...]

💉 In silico serotyping of assemblies.

📥 Inputs:
  database              Database path or keyword (see: `kaptive db list`)
  genomes               Genome assemblies in fasta format; can be compressed

📤 Outputs:
  -o, --out FILE        Write serotyping results as a TSV report to a file (default: stdout)
  -l, --loci [DIR]      Write locus nucleotide fasta files to a directory (default: ./)
  -g, --genes [DIR]     Write gene nucleotide fasta files to a directory (default: ./)
  -p, --proteins [DIR]  Write translation amino-acid fasta files to a directory (default: ./)
  -j, --json [FILE]     Write serialised results to a newline-delimited JSON (default: kaptive_results.jsonl)
  --pha4ge [FILE]       Write PHA4GE-compliant serotyping report to a TSV file (default: kaptive_results.pha4ge)
  --plots [DIR]         Generate interactive locus plots to a directory (default: ./)

🔬 Confidence options:
  --max-other-genes     Typeable if <= other genes (default: 1)
  --min-completeness    Typeable if >= completeness (default: 0.5)
  --below-threshold     Typeable if any genes in locus are below threshold (default: False)

🔧 Other options:
  -t, --threads         Number threads or 0 for all available (default: 0)

🌎 Global options:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress
```
