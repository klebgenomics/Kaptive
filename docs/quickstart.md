# Quickstart

G'day! Welcome to Kaptive! 

If you're looking to get up and running with *in silico* serotyping of genome assemblies as fast as possible, you're in the right place. We'll skip the heavy algorithmic details here and just focus on getting you your results.

## 1. Installation

If you've got Python installed on your machine, you can grab Kaptive straight from PyPI using `pip`:

```bash
pip install kaptive
```

## 2. Download a Database

Kaptive uses beautifully curated databases containing reference surface antigen loci for different species. To see what's available for download, just run:

```bash
kaptive db list
```

To pull down a specific database, use the `install` command with the database name. For example, if you're working with the *Acinetobacter baumannii* K-locus, you'd want `ab_k`:

```bash
kaptive db install ab_k
```

## 3. Serotype Your Genomes!

Once Kaptive and your database are installed, you're ready to roll. Your input genome assemblies should be in FASTA format (and don't worry, gzip-compressed files work perfectly too). 

Let's run Kaptive and save the output to a file called `results.tsv`:

```bash
kaptive type ab_k path/to/genomes/*.fasta -o results.tsv
```

*Pro Tip: You can just use `.` to process all FASTA files in your current directory, or give it explicit paths to your files.*

## 4. Understanding the Output

Kaptive spits out a tab-separated values (TSV) report, which you can easily open up in Excel, Numbers, or any text editor to browse through. 

Here are the most critical columns to keep an eye on in your `results.tsv` file:

* **Assembly**: The name of your input genome file.
* **Match**: The best-matching serotype or locus found in the database (e.g., `KL1`).
* **Confidence**: How confident Kaptive is in this call, ranging from `None` to `Perfect`. A `Perfect`, `Very high`, or `High` confidence generally means you can trust the serotype.
* **Completeness**: The percentage of expected genes for that locus that we actually found in your assembly. A low completeness often means your assembly is a bit fragmented or the locus is truncated.
* **Locus Coordinates**: Exactly where the locus was found within your assembly (Contig, Start, End, and Strand).

For a deeper dive into all the other columns and more advanced outputs, check out the full [Outputs documentation](serotyping/outputs.md). Happy serotyping!
