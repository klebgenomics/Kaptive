
**************************************
Method
**************************************

kaptive assembly
=================
For each input assembly, Kaptive runs the ``kaptive.assembly.typing_pipeline`` which does the following:

#. Aligns locus gene nucleotide sequences to the assembly contig sequences using minimap2.
#. Identifies the best matching locus type using the `Scoring-algorithm`_.
#. Extracts the locus gene sequences from the assembly contig sequences.

.. _Scoring-algorithm:

Scoring algorithm
-------------------
#. For each locus gene, the best alignment is chosen and sorted by locus.
#. For each locus, a chosen alignment metric is summed across genes and weighted to generate a score.
#. The locus with the highest score is chosen as the best match.

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

.. _Locus-reconstruction:

Locus reconstruction
---------------------
After the best matching locus type has been identified, Kaptive will:

#. Align the best matching locus nucleotide sequence to the assembly contig sequences using minimap2.
#. Create pieces of the locus on the contig by merging together ranges of the alignments that are within the distance
   of the largest locus in the database.

.. _Gene-evaluation:

Gene evaluation
---------------------
For each ``GeneResult`` object, Kaptive will:

#. Check whether the gene is **partial** by determining if the gene overlaps the start or end of the contig.
#. Extract the DNA sequence from the assembly contig and translate to amino acid.
#. Perform pairwise alignment to the reference gene amino acid sequence and calculate percent identity.
#. Check for **truncation** by determining if the amino acid sequence length is <95% of the reference gene protein length.
#. Determine whether the gene belongs to the biosynthetic gene cluster (**inside locus**) or not (**outside locus**).

.. note::
 Partial genes are *not* considered for truncation. This prevents false positive truncation calls in
 fragmented assemblies which may otherwise have an impact on phenotype prediction.

.. _Phenotype-prediction:

Phenotype prediction
---------------------
As of Kaptive 3, we have added the ability to predict the resulting phenotype of the assembly. This is similar
to how the *Type* was reported in previous versions, but now includes the ability to predict specific phenotypes
based on known mutations/modifications in a given set of locus genes as defined in the database :ref:`logic file<Phenotype-logic>`.

