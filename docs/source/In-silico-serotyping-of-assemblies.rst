**************************************
*In silico* serotyping of assemblies
**************************************

As with previous versions, typing of bacterial genome assemblies is the primary mode of Kaptive.
Given a Kaptive database and a bacterial genome assembly, this mode  performs 3 main tasks:

* Determines the most likely locus type of the genome assembly.
* Reconstructs the biosynthetic gene cluster from the assembly contig sequences.
* Predicts the corresponding serotype/phenotype of the genome assembly.

.. note::
 As of version 3, Kaptive  no longer supports allelic (*wzi*, *wzc*) typing.

Usage
========

Command-line
--------------

We designed Kaptive 3 to be easier to use on the command-line than previous versions by structuring the program as a
series of sub-commands that follow the general pattern of ``kaptive <mode> <database> <input>``.
To see the full list of commands and options, run ``kaptive --help``.

To perform K-locus typing on a directory of *Klebsiella pneumoniae* assemblies, you would run::

    kaptive assembly kpsc_k assemblies/*.fna

Here we have told Kaptive to perform typing of assemblies with ``assembly`` and used the database keyword
``kpsc_k`` to specify the *Klebsiella pneumoniae* K-locus database.

Other options
^^^^^^^^^^^^^^^^
There are a number of other options for the assembly typing mode:

* ``--min-zscore`` - Minimum zscore for confidence
* ``--min-cov`` - Minimum gene coverage to be used for scoring
* ``--gene-threshold`` - Species-level locus gene identity threshold

We have determined the defaults to be the most appropriate for the majority of use cases, but you can adjust them to
suit your needs.

API
------
Whilst Kaptive isn't designed to be a full API, it is possible to use it as a module in your own Python scripts.
For typing assemblies, you can use the ``kaptive.assembly.typing_pipeline`` function, which takes an assembly path and a
``kaptive.database.Database`` object as input and returns a ``kaptive.typing.TypingResult`` object.

.. code-block:: python

    from kaptive.database import Database, get_database
    from kaptive.assembly import typing_pipeline
    from pathlib import Path

    db = Database.from_genbank(database_path)
    results = [typing_pipeline(assembly, db, threads=8) for assembly in Path('assemblies').glob('*.fna')]

For example, if you wanted to perform K and O locus typing on a single assembly, you could do the following::

    k_db, o_db = get_database('kpsc_k'), get_database('kpsc_o')
    k_db, o_db = Database.from_genbank(k_db), Database.from_genbank(o_db)
    k_results, o_results = typing_pipeline(a, k_db, threads=8), typing_pipeline(a, o_db, threads=8)
    print(k_results.as_table(), o_results.as_table())

Method
========
For each input assembly, Kaptive runs the ``kaptive.assembly.typing_pipeline`` which does the following:

#. Aligns locus gene nucleotide sequences to the assembly contig sequences using minimap2.
#. Identifies the best matching locus type using a novel `scoring algorithm <Scoring-algorithm>`_.
#. Extracts the locus gene sequences from the assembly contig sequences.

Locus reconstruction
---------------------
After the best matching locus type has been identified, Kaptive will attempt to reconstruct the locus biosynthetic gene
cluster from the assembly contig sequences. For each contig where alignments were found, Kaptive will:

#. Sort alignments into whether they belong to the best matching locus (expected) or not (unexpected).
#. Cull all alignments of unexpected genes that overlap with alignments of expected genes.
#. Cull all alignments of unexpected genes that overlap with each other.
#. Create pieces of the locus on the contig by merging together alignment ranges of expected genes that:

   * Are within the distance of the largest locus in the database.
   * Do not extend beyond the alignment ranges of the last expected gene (prevents K-locus matches to the O-locus).
#. ```GeneResult``` objects are created from the remaining alignments then evaluated.

Gene evaluation
^^^^^^^^^^^^^^^^^
For each ```GeneResult``` object, Kaptive will:

#. Check whether the gene is **partial** by determining if the gene overlaps the start or end of the contig.
#. Extract the DNA sequence from the assembly contig and translate to amino acid.
#. Perform pairwise alignment to the reference gene amino acid sequence and calculate percent identity.
#. Check for **truncation** by determining if the amino acid sequence length is >= 95% of the reference gene protein length.
#. Determine whether the gene belongs to the biosynthetic gene cluster (**inside locus**) or not (**outside locus**).

.. note::
 Partial genes are *not* considered for truncation. This prevents false positive truncation calls in
 fragmented assemblies which have an impact on phenotype prediction.

Phenotype prediction
---------------------
As of Kaptive 3, we have added the ability to predict the resulting phenotype of the assembly. This is similar
to how the *Type* was reported in previous versions, but now includes the ability to predict specific phenotypes
based on known mutations/modifications in the locus genes.

Scoring algorithm
-------------------
#. For each locus gene, the best alignment is chosen and sorted by locus.
#. For each locus, a chosen alignment metric is summed across genes and weighted to generate a score.
#. The locus with the highest score is chosen as the best match.
#. The standard deviation is calculated across scores to generate accompanying Z-scores.

The alignment metric can be explicitly set by the flag ``--alignment_metric``; the options are any attribute
of the ``kaptive.alignment.Alignment`` object, such as:

* ``AS`` - alignment score calculated by minimap2 (default)
* ``matching_bases`` - number of matching bases in the alignment
* ``percent_query_coverage`` - the percent of the query sequence covered by the alignment
* ``percent_target_coverage`` - the percent of the target sequence covered by the alignment
* ``percent_identity`` - the percent identity of the alignment

The weighting can be explicitly set by the flag ``--weight_metric``; the options are:

* ``none`` - No weighting
* ``locus_length`` - length of the locus
* ``genes_expected`` - number of genes expected in the locus
* ``genes_found`` - number of genes found in the locus
* ``prop_genes_found`` - number of genes found divided by number of genes expected (default)

Output
==============

Tabular
------------------
The main output of the assembly typing mode is a tab-delimited table of the results with the following columns:

======================================   =====================================================================================================================================
Column name                              Description
======================================   =====================================================================================================================================
Assembly                                 The name of the input assembly, taken from the assembly filename.
Best match locus                         The locus type which most closely matches the assembly.
Best match type                          The predicted serotype/phenotype of the assembly.
Match confidence                         A categorical measure of locus call quality (see confidence reference).
Problems                                 Characters indicating issues with the locus match (see problems reference).
Identity                                 Percent identity of the best matching locus to the assembly.
Coverage                                 Percent coverage of the best matching locus in the assembly.
Length discrepancy                       If the locus was found in a single piece, this is the difference between the locus length and the assembly length.
Expected genes in locus                  A fraction indicating how many of the genes in the best matching locus were found in the locus part of the assembly.
Expected genes in locus, details         Gene names for the expected genes found in the locus part of the assembly.
Missing expected genes                   A string listing the gene names of expected genes that were not found.
Other genes in locus                     The number of unexpected genes (genes from loci other than the best match) which were found in the locus part of the assembly.
Other genes in locus, details            Gene names for the other genes found in the locus part of the assembly.
Expected genes outside locus             A fraction indicating how many of the expected genes which were found in the assembly but not in the locus part of the assembly (usually zero)
Expected genes outside locus, details    Gene names for the expected genes found outside the locus part of the assembly.
Other genes outside locus                The number of unexpected genes (genes from loci other than the best match) which were found outside the locus part of the assembly.
Other genes outside locus, details       Gene names for the other genes found outside the locus part of the assembly.
Truncated genes, details                 Gene names for the truncated genes found in the assembly.
======================================   =====================================================================================================================================

.. note::
 Numbers beside gene names indicate the percent identity and percent coverage of the gene in the assembly.

The default is to print this table to **stdout**.
You can use UNIX redirection operators (``>`` or ``>>``) or the ``-o``/``--out`` flag to write to a file.

If the summary table already exists, Kaptive will append to it (not overwrite it) and suppress the header line.
This allows you to run Kaptive in parallel on many assemblies, all outputting to the same table file.

To disable the tabular output, simply redirect the output to ``/dev/null``.

Fasta
---------
The ``--fasta`` flag produces a fasta file of the region(s) of the assembly which correspond to the best
locus match. This may be a single piece (in cases of a good assembly and a strong match) or it may be in multiple
pieces (in cases of poor assembly and/or a novel locus). The file is named using the output prefix and the assembly name.

The default is to write this file to the current directory with the name: ``{assembly}_kaptive_results.fna``,
however the output directory can be specified after the flag.

JSON
--------
The ``--json`` flag produces a JSON file of the results which allows Kaptive to reconstruct the ``TypingResult`` objects
after a run. Unlike previous version (2 and below), this is a JSON lines file, where each line is a JSON object
representing the results for a single assembly.

The default is to write this file to: ``kaptive_results.json``, however the path can be specified after the flag.
If the file already exists, Kaptive will append to it (not overwrite it).

Diagram
------------------
Kaptive can now produce a visual representation of the locus match in the assembly. This is done using the
``--draw`` flag, which produces a diagram in the format specified by the ``--draw-fmt`` flag (default: png).

The default is to write this file to the current directory with the name: ``{assembly}_kaptive_results.{fmt}``,
however the output directory can be specified after the flag.

Convert
=================

The new ``convert`` command allows you to convert the Kaptive results JSON file into a range of useful formats, including:

* `Tabular <Tabular>`_ output (tsv)
* `Locus nucleotide sequence(s) <Fasta>`_ (fna)
* Locus gene nucleotide sequences (ffn)
* Locus gene amino acid sequences (faa)
* Locus `diagram <Diagram>`_

This means if you didn't want to or forgot to output these files during the initial run, we've got you covered!

Simply run ``kaptive convert <JSON file> <format>`` and the file will be output to the current directory.

