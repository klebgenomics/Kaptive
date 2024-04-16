====================================
How to run
====================================

Quick version (for the impatient)
===================================

Kaptive needs the following input files to run (included in this repository):

* A multi-record Genbank file with your known loci (nucleotide sequences for each whole locus and protein sequences for their genes)

* One or more assemblies in FASTA format

Example command:

.. code-block:: bash

    kaptive.py -a path/to/assemblies/*.fasta -k database.gbk -o output_directory/prefix

For each input assembly file, Kaptive will identify the closest known locus type and report information about the corresponding locus genes.

It generates the following output files:

* A FASTA file for each input assembly with the nucleotide sequences matching the closest locus.
* A table and a JSON file summarising the results for all input assemblies.

Character codes in the output indicate problems with the locus match:

* `?` = the match was not in a single piece, possible due to a poor match or discontiguous assembly.
* `-` = genes expected in the locus were not found.
* `+` = extra genes were found in the locus.
* `*` = one or more expected genes was found but with low identity.

Input files
=============

Assemblies
----------------

Using the ``-a`` (or ``--assembly``) argument, you must provide one or more FASTA files to analyse.
There are no particular requirements about the header formats in these inputs.

Locus references
-------------------

Using the ``-k`` (or ``--k_refs``) argument, you must provide a Genbank file containing one record for each known locus.

This input Genbank has the following requirements:

* The ``source`` feature must contain a ``note`` qualifier which begins with a label such as ``K locus:``.
  Whatever follows is used as the locus name reported in the Kaptive output. The label is automatically determined,
  and any consistent label ending in a colon will work. However, the user can specify exactly which label to use with
  ``--locus_label``, if desired.

* The ``source`` feature may optionally contain a ``note`` qualifier which begins with a label such as
  ``K type:`` that specifies the serotype (phenotype) associated with the locus (is known). In cases where only some loci
  are associated with known serotypes we recommend adding a ``note`` such as ``K type: unknown``. If no ``type`` notes
  are specified for any loci, the Kaptive will list them as ``unknown`` in the output. (Kaptive v2.0+)

* Any locus gene should be annotated as ``CDS`` features. All ``CDS`` features will be used and any other type of
  feature will be ignored.

* If the gene has a name, it should be specified in a ``gene`` qualifier. This is not required, but if absent the gene
  will only be named using its numbered position in the locus.

Example piece of input Genbank file:

.. code-block::

    source          1..23877
                    /organism="Klebsiella pneumoniae"
                    /mol_type="genomic DNA"
                    /note="K locus: KL1"
                    /note="K type: K1"
    CDS             1..897
                    /gene="galF"



Allelic typing
------------------

You can also supply Kaptive with a FASTA file of gene alleles using the ``-g`` (or ``--allelic-typing``) argument.
For example, ``wzi_wzc_db.fasta`` (included with Kaptive) contains *Klebsiella* wzi and wzc alleles. This file must be
formatted as an `SRST2 <https://github.com/katholt/srst2>`_ database with integers for allele names. If used, Kaptive
will report the number of the best allele for each type gene within the found locus. If there is no perfect match,
Kaptive reports the best match and adds a ``*`` to the allele number. To be clear, Kaptive will only search for the genes
in the locus sequence, not in the entire genome (check out `this tool <https://github.com/rrwick/SRST2-table-from-assemblies>`_
if an entire-genome search is needed).


Standard output
=================

Basic
--------

Kaptive will write a simple line to stdout for each assembly:

* the assembly name
* the best locus match
* character codes for any match problems

Example (no type genes supplied):

.. code-block::

    assembly_1: K2*
    assembly_2: K4
    assembly_3: KL17?-*

Example (with type genes):

.. code-block::

    assembly_1: K2*, wzc=2, wzi=2
    assembly_2: K4, wzc=1, wzi=127
    assembly_3: KL17?-*, wzc=18*, wzi=137*


Verbose
----------

If run without the ``-v`` or ``--verbose`` option, Kaptive will give detailed information about each assembly including:

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


Output files
==============

Summary table
------------------

Kaptive produces a single tab-delimited table summarising the results of all input assemblies. It has the following columns:

* **Assembly**: the name of the input assembly, taken from the assembly filename.
* **Best match locus**: the locus type which most closely matches the assembly, based on BLAST coverage.
* **Best match type**: the predicted serotype associated with the best match locus, as specified in the reference database (Kaptive v2.0+)
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

To disable the table file output, run Kaptive with ``--no_table``.


JSON
========

Kaptive also outputs its results in a JSON file which contains all information from the above table, as well as more detail about BLAST results and reference sequences.

To disable JSON output, run Kaptive with ``--no_json``.


Locus matching sequences
============================

For each input assembly, Kaptive produces a fasta file of the region(s) of the assembly which correspond to the best locus match. This may be a single piece (in cases of a good assembly and a strong match) or it may be in multiple pieces (in cases of poor assembly and/or a novel locus). The file is named using the output prefix and the assembly name.

To disable these output files, run Kaptive with ``--no-seq-out``.

Advanced options
==================

Each of these options has a default and is not required on the command line, but can be adjusted if desired:

* ``--start-end-margin``: Kaptive tries to identify whether the start and end of a locus are present in an assembly and in the same contig. This option allows for a bit of wiggle room in this determination. For example, if this value is 10 (the default), a locus match that is missing the first 8 base pairs will still count as capturing the start of the locus. If set to zero, then the BLAST hit(s) must extend to the very start and end of the locus for Kaptive to consider the match complete.
* ``--min-gene-cov``: the minimum required percent coverage for the gene BLAST search via tBLASTn. For example if this value is 90 (the default), then a gene BLAST hit which only covers 85% of the gene will be ignored. Using a lower value will allow smaller pieces of genes to be included in the results.
* ``--min-gene-id``: the mimimum required percent identity for the gene BLAST search via tBLASTn. For example if this value is 80 (the default), then a gene BLAST hit which has only 65% amino acid identity will be ignored. A lower value will allow for more distant gene hits to be included in the results (possibly resulting in more genes in the 'Other genes outside locus' category). A higher value will make Kaptive only accept very close gene hits (possibly resulting in low-identity locus genes not being found and included in 'Missing expected genes').
* ``--low-gene-id``: the percent identity threshold for what counts as a low identity match in the gene BLAST search. This only affects whether or not the `*` character is included in the 'Problems'. Default is 95.
* ``--min-assembly-piece``: the smallest piece of the assembly (measured in bases) that will be included in the output FASTA files. For example, if this value is 100 (the default), then a 50 bp match between the assembly and the best matching locus reference will be ignored.
* ``--gap-fill-size``: the size of assembly gaps to be filled in when producing the output FASTA files. For example, if this value is 100 (the default) and an assembly has two separate locus BLAST hits which are only 50 bp apart in a contig, they will be merged together into one sequence for the output FASTA. But if the two BLAST hits were 150 bp apart, they will be included in the output FASTA as two separate sequences. A lower value will possibly result in more fragmented output FASTA sequences. A higher value will possibly result in more sequences being included in the locus output.


SLURM jobs
============

If you are running this script on a cluster using `SLURM <http://slurm.schedmd.com/>`_, then you can make use of the
extra script: ``kaptive_slurm.py``. This will create one SLURM job for each assembly so the jobs can run in parallel.
All simultaneous jobs can write to the same output table. It may be necessary to modify this script to suit the details
of your cluster.
