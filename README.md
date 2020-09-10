<p align="center"><img src="extras/kaptive_logo.png" alt="Kaptive" width="400"></p>


Kaptive reports information about surface polysaccharide loci for _Klebsiella_ and _Acinetobacter baumannii_ genome assemblies. You can also run a graphical version of Kaptive via [this web interface](http://kaptive.holtlab.net/) ([source code](https://github.com/kelwyres/Kaptive-Web)).

Given a novel genome and a database of known loci (K, O or OC), Kaptive will help a user to decide whether their sample has a known or novel locus. It carries out the following for each input assembly:
* BLAST for all known locus nucleotide sequences (using `blastn`) to identify the best match ('best' defined as having the highest coverage).
* Extract the region(s) of the assembly which correspond to the BLAST hits (i.e. the locus sequence in the assembly) and save it to a FASTA file.
* BLAST for all known locus genes (using `tblastn`) to identify which expected genes (genes in the best matching locus) are present/missing and whether any unexpected genes (genes from other loci) are present.
* Output a summary to a table file.

In cases where your input assembly closely matches a known locus, Kaptive should make that obvious. When your assembly has a novel type, that too should be clear. However, Kaptive cannot reliably extract or annotate locus sequences for totally novel types – if it indicates a novel locus is present then extracting and annotating the sequence is up to you! Very poor assemblies can confound the results, so be sure to closely examine any case where the locus sequence in your assembly is broken into multiple pieces.
If you think you have found a novel locus that should be added to one of the databases distributed with Kaptive please [contact us](mailto:kaptive.typing@gmail.com).

For citation info and details about Kaptive, Kaptive Web and the locus databases, see [our papers](#citation) below.


## Table of Contents

* [Quick version (for the impatient)](#quick-version-for-the-impatient)
* [Installation](#installation)
* [Input files](#input-files)
* [Standard output](#standard-output)
    * [Basic](#basic)
    * [Verbose](#verbose)
* [Output files](#output-files)
    * [Summary table](#summary-table)
    * [JSON](#json)
    * [Locus matching sequences](#locus-matching-sequences)
* [Example results and interpretation](#example-results-and-interpretation)
    * [Very close match](#very-close-match)
    * [More distant match](#more-distant-match)
    * [Broken assembly](#broken-assembly)
    * [Poor match](#poor-match)
* [Advanced options](#advanced-options)
* [SLURM jobs](#slurm-jobs)
* [Databases distributed with Kaptive](#databases-distributed-with-kaptive)
    * [_Klebsiella_ K locus databases](#klebsiella-k-locus-databases)
    * [_Klebsiella_ O locus database](#klebsiella-o-locus-database)
    * [_Acinetobacter baumannii_ K and OC locus databases](#acinetobacter-baumanii-k-and-oc-locus-databases)
* [FAQs](#faqs)
* [Citation](#citation)
* [License](#license)


## Quick version (for the impatient)

Kaptive needs the following input files to run (included in this repository):
* A multi-record Genbank file with your known loci (nucleotide sequences for each whole locus and protein sequences for their genes)
* One or more assemblies in FASTA format

Example command:

`kaptive.py -a path/to/assemblies/*.fasta -k database.gbk -o output_directory/prefix`

For each input assembly file, Kaptive will identify the closest known locus type and report information about the corresponding locus genes.

It generates the following output files:
* A FASTA file for each input assembly with the nucleotide sequences matching the closest locus
* A table summarising the results for all input assemblies

Character codes in the output indicate problems with the locus match:
* `?` = the match was not in a single piece, possible due to a poor match or discontiguous assembly.
* `-` = genes expected in the locus were not found.
* `+` = extra genes were found in the locus.
* `*` = one or more expected genes was found but with low identity.


## Installation

Kaptive should work on both Python 2 and 3, but we run/test it on Python 3 and recommend you do the same.


#### Clone and run

Kaptive is a single Python script, so you can simply clone (or download) from GitHub and run it. Kaptive depends on [Biopython](http://biopython.org/wiki/Main_Page), so make sure it's installed (either `pip3 install biopython` or read detailed instructions [here](http://biopython.org/DIST/docs/install/Installation.html)).

```
git clone https://github.com/katholt/Kaptive
Kaptive/kaptive.py -h
```

#### Install with pip

Alternatively, you can install Kaptive using [pip](https://pip.pypa.io/en/stable/). This will take care of the Biopython requirement (if necessary) and put the `kaptive.py` script in your PATH for easy access. Pip installing will _not_ provide the reference databases, so you'll need to download them separately from [here](https://github.com/katholt/Kaptive/tree/master/reference_database).

```
pip3 install --user git+https://github.com/katholt/Kaptive
kaptive.py -h
```

#### Other dependencies

Regardless of how you download/install Kaptive, it requires that [BLAST+](http://www.ncbi.nlm.nih.gov/books/NBK279690/) is available on the command line (specifically the commands `makeblastdb`, `blastn` and `tblastn`). BLAST+ can usually be easily installed using a package manager such as [Homebrew](http://brew.sh/) (on Mac) or [apt-get](https://help.ubuntu.com/community/AptGet/Howto) (on Ubuntu and related Linux distributions). Some later versions of BLAST+ have been associated with sporadic crashes when running tblastn with multiple threads; to avoid this problem we recommend running Kaptive with BLAST+ v 2.3.0 or using the "--threads 1" option (see below for full command argument details).


## Input files

#### Assemblies

Using the `-a` (or `--assembly`) argument, you must provide one or more FASTA files to analyse. There are no particular requirements about the header formats in these inputs.

#### Locus references

Using the `-k` (or `--k_refs`) argument, you must provide a Genbank file containing one record for each known locus.

This input Genbank has the following requirements:
* The `source` feature must contain a `note` qualifier which begins with a label such as 'K locus:'. Whatever follows is used as the locus name. The label is automatically determined, and any consistent label ending in a colon will work. However, the user can specify exactly which label to use with `--locus_label`, if desired.
* Any locus gene should be annotated as `CDS` features. All `CDS` features will be used and any other type of feature will be ignored.
* If the gene has a name, it should be specified in a `gene` qualifier. This is not required, but if absent the gene will only be named using its numbered position in the locus.

Example piece of input Genbank file:
```
source          1..23877
                /organism="Klebsiella pneumoniae"
                /mol_type="genomic DNA"
                /note="K locus: K1"
CDS             1..897
                /gene="galF"
```

#### Allelic typing

You can also supply Kaptive with a FASTA file of gene alleles using the `-g` (or `--allelic_typing`) 
argument. For example, `wzi_wzc_db.fasta` (included with Kaptive) contains wzi and wzc alleles. This file must be formatted as an [SRST2](https://github.com/katholt/srst2) database with integers for allele names.

If used, Kaptive will report the number of the best allele for each type gene. If there is no perfect match, Kaptive reports the best match and adds a `*` to the allele number.


## Standard output

#### Basic

Kaptive will write a simple line to stdout for each assembly:
* the assembly name
* the best locus match
* character codes for any match problems

Example (no type genes supplied):
```
assembly_1: K2*
assembly_2: K4
assembly_3: KL17?-*
```

Example (with type genes):
```
assembly_1: K2*, wzc=2, wzi=2
assembly_2: K4, wzc=1, wzi=127
assembly_3: KL17?-*, wzc=18*, wzi=137*
```

#### Verbose

If run without the `-v` or `--verbose` option, Kaptive will give detailed information about each assembly including:
* Which locus reference best matched the assembly
* Information about the nucleotide sequence match between the assembly and the best locus reference:
  * % Coverage and % identity
  * Length discrepancy (only available if assembled locus match is in one piece)
  * Contig names and coordinates for matching sequences
* Details about found genes:
  * Whether they were expected or unexpected
  * Whether they were found inside or outside the locus matching sequence
  * % Coverage and % identity
  * Contig names and coordinates for matching sequences
* Best alleles for each type gene (if the user supplied a type gene database)


## Output files

#### Summary table

Kaptive produces a single tab-delimited table summarising the results of all input assemblies. It has the following columns:
* **Assembly**: the name of the input assembly, taken from the assembly filename.
* **Best match locus**: the locus type which most closely matches the assembly, based on BLAST coverage.
* **Match confidence**: a categorical measure of match quality:
  * `Perfect` = the locus was found in a single piece with 100% coverage and 100% identity.
  * `Very high` = the locus was found in a single piece with ≥99% coverage and ≥95% identity, with no missing genes and no extra genes.
  * `High` = the locus was found in a single piece with ≥99% coverage, with ≤ 3 missing genes and no extra genes.
  * `Good` = the locus was found in a single piece or with ≥95% coverage, with ≤ 3 missing genes and ≤ 1 extra genes.
  * `Low` = the locus was found in a single piece or with ≥90% coverage, with ≤ 3 missing genes and ≤ 2 extra genes.
  * `None` = did not qualify for any of the above.
* **Problems**: characters indicating issues with the locus match. An absence of any such characters indicates a very good match.
  * `?` = the match was not in a single piece, possible due to a poor match or discontiguous assembly.
  * `-` = genes expected in the locus were not found.
  * `+` = extra genes were found in the locus.
  * `*` = one or more expected genes was found but with low identity.
* **Coverage**: the percent of the locus reference which BLAST found in the assembly.
* **Identity**: the nucleotide identity of the BLAST hits between locus reference and assembly.
* **Length discrepancy**: the difference in length between the locus match and the corresponding part of the assembly. Only available if the locus was found in a single piece (i.e. the `?` problem character is not used).
* **Expected genes in locus**: a fraction indicating how many of the genes in the best matching locus were found in the locus part of the assembly.
* **Expected genes in locus, details**: gene names and percent identity (from the BLAST hits) for the expected genes found in the locus part of the assembly.
* **Missing expected genes**: a string listing the gene names of expected genes that were not found.
* **Other genes in locus**: the number of unexpected genes (genes from loci other than the best match) which were found in the locus part of the assembly.
* **Other genes in locus, details**: gene names and percent identity (from the BLAST hits) for the other genes found in the locus part of the assembly.
* **Expected genes outside locus**: the number of expected genes which were found in the assembly but not in the locus part of the assembly (usually zero)
* **Expected genes outside locus, details**: gene names and percent identity (from the BLAST hits) for the expected genes found outside the locus part of the assembly.
* **Other genes outside locus**: the number of unexpected genes (genes from loci other than the best match) which were found outside the locus part of the assembly.
* **Other genes outside locus, details**: gene names and percent identity (from the BLAST hits) for the other genes found outside the locus part of the assembly.
* One column for each type gene (if the user supplied a type gene database)

If the summary table already exists, Kaptive will append to it (not overwrite it). This allows you to run Kaptive in parallel on many assemblies, all outputting to the same table file.

To disable the table file output, run Kaptive with `--no_table`.


#### JSON

Kaptive also outputs its results in a JSON file which contains all information from the above table, as well as more detail about BLAST results and reference sequences.

To disable JSON output, run Kaptive with `--no_json`.


#### Locus matching sequences

For each input assembly, Kaptive produces a Genbank file of the region(s) of the assembly which correspond to the best locus match. This may be a single piece (in cases of a good assembly and a strong match) or it may be in multiple pieces (in cases of poor assembly and/or a novel locus). The file is named using the output prefix and the assembly name.

To these output files, run Kaptive with `--no_seq_out`.


## Example results and interpretation

These examples show what Kaptive's results might look like in the output table. The gene details columns of the table have been excluded for brevity, as they can be quite long.

#### Very close match

Assembly | Best match locus | Problems | Coverage | Identity | Length discrepancy | Expected genes in locus | Expected genes in locus, details | Missing expected genes | Other genes in locus | Other genes in locus, details | Expected genes outside locus | Expected genes outside locus, details | Other genes outside locus | Other genes outside locus, details
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
assembly_1 | K1 |  | 99.94% | 99.81% | -22 bp | 20 / 20 (100%) | ... |  | 0 |  | 0 |  | 2 | ...

This is a case where our assembly very closely matches a known K locus type. There are no characters in the 'Problems' column, the coverage and identity are both high, the length discrepency is low, and all expected genes were found with high identity. A couple of other low-identity K locus genes hits were elsewhere in the assembly, but that's not abnormal and no cause for concern.

Overall, this is a nice, solid match for K1.

#### More distant match

Assembly | Best match locus | Problems | Coverage | Identity | Length discrepancy | Expected genes in locus | Expected genes in locus, details | Missing expected genes | Other genes in locus | Other genes in locus, details | Expected genes outside locus | Expected genes outside locus, details | Other genes outside locus | Other genes outside locus, details
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
assembly_2 | K1 | * | 99.84% | 95.32% | +97 bp | 20 / 20 (100%) | ... |  | 0 |  | 0 |  | 2 | ...

This case shows an assembly that also matches the K1 locus sequence, but not as closely as our previous case. The `*` character indicates that one or more of the expected genes falls below the identity threshold (default 95%). The 'Expected genes in K locus, details' columns, excluded here for brevity, would show the identity for each gene.

Our sample still almost certainly has a K locus type of K1, but it has diverged a bit more from our K1 reference, possibly due to mutation and/or recombination.

#### Broken assembly

Assembly | Best match locus | Problems | Coverage | Identity | Length discrepancy | Expected genes in locus | Expected genes in locus, details | Missing expected genes | Other genes in locus | Other genes in locus, details | Expected genes outside locus | Expected genes outside locus, details | Other genes outside locus | Other genes outside locus, details
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
assembly_3 | K2 | ?- | 99.95% | 98.38% | n/a | 17 / 18 (94.4%) | ... | K2-CDS17-manB | 0 |  | 0 |  | 1 | ...

Here is a case where our assembly matched a known K locus type well (high coverage and identity) but with a couple of problems. First, the `?` character indicates that the K locus sequence was not found in one piece in the assembly. Second, one of the expected genes (K2-CDS17-manB) was not found in the gene BLAST search.

In cases like this, it is worth examining the case in more detail outside of Kaptive. For this example, such an examination revealed that the assembly was poor (broken into many small pieces) and the _manB_ gene happened to be split between two contigs. So the _manB_ gene isn't really missing, it's just broken in two. Our sample most likely is a very good match for K2, but the poor assembly quality made it difficult for Kaptive to determine that automatically.

#### Poor match

Assembly | Best match locus | Problems | Coverage | Identity | Length discrepancy | Expected genes in locus | Expected genes in locus, details | Missing expected genes | Other genes in locus | Other genes in locus, details | Expected genes outside locus | Expected genes outside locus, details | Other genes outside locus | Other genes outside locus, details
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
assembly_4 | K3 | ?-* | 77.94% | 83.60% | n/a | 15 / 20 (75%) | ... | ... | 0 |  | 0 |  | 5 | ...

In this case, Kaptive did not find a close match to any known K locus sequence. The best match was to K3, but BLAST only found alignments for 78% of the K3 sequence, and only at 84% nucleotide identity. Five of the twenty K3 genes were not found, and the 15 which were found had low identity. The assembly sequences matching K3 did not come in one piece (indicated by `?`), possibly due to assembly problems, but more likely due to the fact that our sample is not in fact K3 but rather has some novel K locus that was not in our reference inputs.

A case such as this demands a closer examination outside of Kaptive. It is likely a novel K locus type, and you may wish to extract and annotate the K locus sequence from the assembly.


## Advanced options

Each of these options has a default and is not required on the command line, but can be adjusted if desired:

* `--start_end_margin`: Kaptive tries to identify whether the start and end of a locus are present in an assembly and in the same contig. This option allows for a bit of wiggle room in this determination. For example, if this value is 10 (the default), a locus match that is missing the first 8 base pairs will still count as capturing the start of the locus. If set to zero, then the BLAST hit(s) must extend to the very start and end of the locus for Kaptive to consider the match complete.
* `--min_gene_cov`: the minimum required percent coverage for the gene BLAST search. For example if this value is 90 (the default), then a gene BLAST hit which only covers 85% of the gene will be ignored. Using a lower value will allow smaller pieces of genes to be included in the results.
* `--min_gene_id`: the mimimum required percent identity for the gene BLAST search. For example if this value is 80 (the default), then a gene BLAST hit which has only 65% amino acid identity will be ignored. A lower value will allow for more distant gene hits to be included in the results (possibly resulting in more genes in the 'Other genes outside locus' category). A higher value will make Kaptive only accept very close gene hits (possibly resulting in low-identity locus genes not being found and included in 'Missing expected genes').
* `--low_gene_id`: the percent identity threshold for what counts as a low identity match in the gene BLAST search. This only affects whether or not the `*` character is included in the 'Problems'. Default is 95.
* `--min_assembly_piece`: the smallest piece of the assembly (measured in bases) that will be included in the output FASTA files. For example, if this value is 100 (the default), then a 50 bp match between the assembly and the best matching locus reference will be ignored.
* `--gap_fill_size`: the size of assembly gaps to be filled in when producing the output FASTA files. For example, if this value is 100 (the default) and an assembly has two separate locus BLAST hits which are only 50 bp apart in a contig, they will be merged together into one sequence for the output FASTA. But if the two BLAST hits were 150 bp apart, they will be included in the output FASTA as two separate sequences. A lower value will possibly result in more fragmented output FASTA sequences. A higher value will possibly result in more sequences being included in the locus output.


## SLURM jobs

If you are running this script on a cluster using [SLURM](http://slurm.schedmd.com/), then you can make use of the extra script: `kaptive_slurm.py`. This will create one SLURM job for each assembly so the jobs can run in parallel. All simultaneous jobs can write to the same output table. It may be necessary to modify this script to suit the details of your cluster.



## Databases distributed with Kaptive

The Kaptive repository contains _Klebsiella_ K-locus and O-locus databases plus _A. baumannii_ K-locus and OC-locus databases in the [reference_database](https://github.com/katholt/Kaptive/tree/master/reference_database) directory, but you can run Kaptive with any appropriately formatted database of your own.

The databases were developed and curated by [Kelly Wyres](https://holtlab.net/kelly-wyres/) (_Klebsiella_) and [Johanna Kenyon](https://research.qut.edu.au/infectionandimmunity/projects/bacterial-polysaccharide-research/) (_A. baumannii_).

If you have a locus database that you would like to be added to Kaptive for use by yourself and others in the community, please get in touch via the [issues page](https://github.com/katholt/Kaptive/issues) or [email](mailto:kaptive.typing@gmail.com) . Similarly, if you have identified new locus variants not currently in the existing databases, let us know!

#### _Klebsiella_ K locus databases

The _Klebsiella_ K locus primary reference database (`Klebsiella_k_locus_primary_reference.gbk`) comprises full-length (_galF_ to _ugd_) annotated sequences for each distinct _Klebsiella_ K locus, where available:
* KL1 - KL77 correspond to the loci associated with each of the 77 serologically defined K-type references.
* KL101 and above are defined from DNA sequence data on the basis of gene content.

Note that insertion sequences (IS) are excluded from this database since we assume that the ancestral sequence was likely IS-free and IS transposase genes are not specific to the K locus.
Synthetic IS-free K locus sequences were generated for K loci for which no naturally occurring IS-free variants have been identified to date.

The variants database (`Klebsiella_k_locus_variant_reference.gbk`) comprises full-length annotated sequences for variants of the distinct loci:
* IS variants are named as KLN -1, -2 etc e.g. KL15-1 is an IS variant of KL15.
* Deletion variants are named KLN-D1, -D2 etc e.g. KL15-D1 is a deletion variant of KL15.
Note that KL156-D1 is included in the primary reference database since no full-length version of this locus has been identified to date. 

We recommend screening your data with the primary reference database first to find the best-matching K locus type. If you have poor matches or are particularly interested in detecting variant loci you should try the variant database.  
WARNING: If you use the variant database please inspect your results carefully and decide for yourself what constitutes a confident match! Kaptive is not optimised for accurate variant detection. 

Database versions:
* Kaptive releases v0.5.1 and below include the original _Klebsiella_ K locus databases, as described in [Wyres, K. et al. Microbial Genomics (2016).](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102)
* Kaptive v0.6.0 and above include four novel primary _Klebsiella_ K locus references defined on the basis of gene content (KL162-KL165) in this [paper.](https://www.biorxiv.org/content/10.1101/557785v1)
* Kaptive v0.7.1 and above contain updated versions of the KL53 and KL126 loci (see table below for details). The updated KL126 locus sequence will be described in McDougall, F. et al. 2020. _Klebsiella pneumoniae_ diversity and detection of _Klebsiella africana_ in Australian Fruit Bats (_Pteropus policephalus_). _In prep._
* Kaptive v0.7.2 and above include a novel primary _Klebsiella_ K locus reference defined on the basis of gene content (KL166), which will be described in Li, M. et al. 2020. Characterization of clinically isolated hypermucoviscous _Klebsiella pneumoniae_ in Japan. _In prep._
* Kaptive v0.7.3 and above include four novel primary _Klebsiella_ K locus references defined on the basis of gene content (KL167-KL170), which will be described in Gorrie, C. et al. 2020. Opportunity and diversity: A year of _Klebsiella pneumoniae_ infections in hospital. _In prep._


Changes to the _Klebsiella_ K locus primary reference database:

| Locus  | Change | Reason | Date of change | Kaptive version no. |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| KL53  | Annotation update: _wcaJ_ changed to _wbaP_ | Error in original annotation | 21 July 2020 | v 0.7.1 | 
| KL126  | Sequence update: new sequence from isolate FF923 includes _rmlBADC_ genes between _gnd_ and _ugd_ | Assembly scaffolding error in original sequence from isolate A-003-I-a-1 | 21 July 2020 | v 0.7.1 |

#### _Klebsiella_ O locus database

The _Klebsiella_ O locus database (`Klebsiella_o_locus_primary_reference.gbk`) contains annotated sequences for 12 distinct _Klebsiella_ O loci.

O locus classification requires some special logic, as the O1 and O2 serotypes contain the same locus genes. It is two additional genes elsewhere in the chromosome (_wbbY_ and _wbbZ_) which results in the O1 antigen. Kaptive therefore looks for these genes to properly call an assembly as either O1 or O2. When only one of the two additional genes can be found, the result is ambiguous and Kaptive will report a locus type of O1/O2.

Read more about the O locus and its classification here: [The diversity of _Klebsiella_ pneumoniae surface polysaccharides](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320592/).

Database versions:
* Kaptive v0.4.0 and above include the original version of the _Klebsiella_ O locus database, as described in [Wick, R. et al. J Clin Microbiol (2019).](http://jcm.asm.org/content/56/6/e00197-18) 


#### _Acinetobacter baunannii_ K and OC locus databases

The _A. baumannii_ K (capsule) locus reference database (`Acinetobacter_baumannii_k_locus_primary_reference.gbk`) contains annotated sequences for 92 distinct K loci.
The _A. baumannii_ OC (lipooligosaccharide outer core) locus reference database (`Acinetobacter_baumannii_OC_locus_primary_reference.gbk`) contains annotated sequences for 12 distinct OC loci.

WARNING: These databases have been developed and tested specifically for _A. baumannii_ and may not be suitable for screening other _Acinetobacter_ species. You can check that your assembly is a true _A. baumannii_ by screening for the _oxaAB_ gene e.g. using blastn.

 Database versions:
* Kaptive v0.7.0 and above include the original _A. baumannii_ K and OC locus databases, as described in [Wyres, KL. et al. Microbial Genomics, 2020.](https://doi.org/10.1099/mgen.0.000339)



## FAQs

#### Why are there locus genes found outside the locus?

For _Klebsiella_ K loci in particular, a number of the K-locus genes are orthologous to genes outside of the K-locus region of the genome. E.g the _Klebsiella_ K-locus <i>man</i> and <i>rml</i> genes have orthologues in the LPS (lipopolysacharide) locus; so it is not unusual to find a small number of genes "outside" the locus.
However, if you have a large number of genes (>5) outside the locus it may mean that there is a problem with the locus match, or that your assembly is very fragmented or contaminated (contains more than one sample).

#### How can my sample be missing locus genes when it has a full-length, high identity locus match?

Kaptive uses 'tblastn' to screen for the presence of each locus gene with a coverage threshold of 90%. A single non-sense mutation or small indel in the centre of a gene will interrupt the 'tblastn' match and cause it to fall below the 90% threshold. However, such a small change has only a minor effect on the nucleotide 'blast' match across the full locus.

#### Why does the _Klebsiella_ K-locus region of my sample contain a <i>ugd</i> gene matching another locus?

A small number of the original _Klebsiella_ K locus references are truncated, containing only a partial <i>ugd</i> sequence. The reference annotations for these loci do not include <i>ugd</i>, so are not identified by the 'tblastn' search. Instead <b>Kaptive</b> reports the closest match to the partial sequence (if it exceeds the 90% coverage threshold). 

#### Why has the best matching locus changed after I reran my analysis with an updated version of the database? ####

The databases are updated as novel loci are discovered and curated. If your previous match had a confidence call of 'Low' or 'None' but your new match has higher confidence, this indicates that your genome contains a locus that was absent in the older version of the database! So nothing to worry about here.

But what if your old match and your new match have 'Good' or better confidence levels?

If your old match had 'Perfect' or 'Very High' confidence, please post an issue to the issues page, as this may indicate a problem with the new database!

If your old match had 'Good' or 'High' confidence please read on...

Polysaccharide loci are subject to frequent recombinations and rearrangements, which generates new variants. As a result, a small number of pairs of loci share large regions of homology e.g. the _Klebsiella_ K-locus KL170 is very similar to KL101, and in fact seems to be a hybrid of KL101 plus a small region from KL106. 
Kaptive can accurately distinguish the KL101 and KL170 loci when it is working with high quality genome assemblies, but this task is much trickier if the assembly is fragmented. This means that matches to KL101 that were reported using an early version of the K-locus database might be reported as KL170 when using a later version of the database.
However, this should only occur in instances where the K-locus is fragmented in the genome assembly and in that case Kaptive will have indicated 'problems' with the matches (e.g. '?' indicating fragmented assembly or '-' indicating that an expected gene is missing), and the corresponding confidence level will be at the lower end of the scale (i.e. 'Good' or 'High', but not 'Very High' or 'Perfect').
You may want to try to figure out the correct locus manually, e.g. using [Bandage](https://rrwick.github.io/Bandage/) to BLAST the corresponding loci in your genome assembly graph. 


## Citation

If you use Kaptive and/or the _Klebsiella_ K locus database in your research, please cite this paper:
[Identification of _Klebsiella_ capsule synthesis loci from whole genome data. Microbial Genomics (2016).](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102)

If you use [Kaptive Web](http://kaptive.holtlab.net/) and/or the _Klebsiella_ O locus database in your research, please cite this paper:
[Kaptive Web: user-friendly capsule and lipopolysaccharide serotype prediction for _Klebsiella_ genomes. Journal of Clinical Microbiology (2018).](http://jcm.asm.org/content/56/6/e00197-18)

If you use the _A. baumannii_ K or OC locus database(s) in your research please cite this paper:
[Identification of _Acinetobacter baumannii_ loci for capsular polysaccharide (KL) and lipooligosaccharide outer core (OCL) synthesis in genome assemblies using curated reference databases compatible with Kaptive. Microbial Genomics (2020).](https://doi.org/10.1099/mgen.0.000339)  
Lists of papers describing each of the individual _A. baumannii_ reference loci can be found [here](https://github.com/katholt/Kaptive/tree/master/extras).


## License

GNU General Public License, version 3

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.495263.svg)](https://doi.org/10.5281/zenodo.495263)
