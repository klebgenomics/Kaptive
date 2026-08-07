---
title: type (assembly)
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/syringe
categories:
  - Serotyping
---

💉 In silico serotyping of genome assemblies.

Aliases:
    assembly


```text
[1;36musage: [0m[1;36m[0mkaptive type [1;36m[options][0m database genomes [genomes ...]

[1m💉 In silico serotyping of genome assemblies.

Aliases:
    assembly
[0m

[1;36m[1m📥 Inputs[0m[0m:
  database              Database path or keyword (see: `kaptive db list`)
  genomes               Genome assemblies in fasta format; can be compressed

[1;36m[1m📤 Outputs[0m[0m:
  -o, --out FILE        Write serotyping results as a TSV report to a file (default: stdout)
  -l, --loci [DIR]      Write locus nucleotide fasta files to a directory (default: ./)
  -g, --genes [DIR]     Write gene nucleotide fasta files to a directory (default: ./)
  -p, --proteins [DIR]  Write translation amino-acid fasta files to a directory (default: ./)
  -j, --json [FILE]     Write serialised results to a newline-delimited JSON (default: kaptive_results.jsonl)
  --pha4ge [FILE]       Write PHA4GE-compliant serotyping report to a TSV file (default: kaptive_results.pha4ge)
  --plots [DIR]         Generate interactive locus plots to a directory (default: ./)

[1;36m[1m🔬 Confidence options[0m[0m:
  --max-other-genes     Typeable if <= other genes (default: 1)
  --min-completeness    Typeable if >= completeness (default: 0.5)
  --below-threshold     Typeable if any genes in locus are below threshold (default: False)

[1;36m[1m🔧 Other options[0m[0m:
  -t, --threads         Number threads or 0 for all available (default: 0)
  --partial-edge-tolerance 
                        Tolerance in bases from contig edge to call a partial gene (default: 5)

[1;36m[1m🌎 Global options[0m[0m:
  -h, --help            show this help message and exit
  -V, --verbose         Enable verbose output/progress
```
