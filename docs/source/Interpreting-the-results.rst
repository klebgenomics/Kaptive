**************************************
Interpreting the results
**************************************

The four most important columns in the Kaptive :ref:`tabular <Tabular>` output:

* ``Best match locus`` indicating the locus that is best supported by the available sequence data
* ``Best match type`` indicating the predicted phenotype based on that associated with the ``Best match locus`` and taking into account any special phenotype logic e.g. where a known essential gene is truncated, Kaptive will report the ``Best match type`` as 'Capsule null'. 
* ``Match confidence`` indicates a qualitative level of confidence that the reported ``Best match locus`` is correct (see below).
* ``Problems`` indicates the features of the assembly locus that may impact ``Match confidence``

.. note::
  The ``Match confidence`` relates to the locus match, and is not a direct indication of confidence in the ``Best match type``.  


.. _Confidence-score:

Confidence score  
=================
Kaptive will indicate the best matching locus and its confidence in the locus match.


Typeable loci
---------------------

The locus was found in a single piece in the query assembly with no genes below the minimum translated identity
according to the :ref:`locus thresholds <Locus definition>` and:

* no missing genes (as determined by at least 50% nucleotide mapping coverage)
* no unexpected genes (genes from other loci) inside the locus region of the assembly

**OR**

The locus was found in more than one piece in the query assembly with no genes below the minimum translated identity
according to the :ref:`locus thresholds <Locus definition>` and:

* no more than 1 missing gene  
* no more than 1 unexpected gene (genes from other loci) inside the locus region of the assembly

These criteria were designed in consideration of the locus definition rules (i.e. that each locus represents a unique set of genes defined at a given minimum translated identity threshold) and following systematic analysis of Kaptive outputs for draft genome assemblies compared against manually confirmed loci determined from matched completed genomes.

We allow some flexibility with regards to missing genes or additional genes found within the locus when this region of the query assembly is fragmented, because it can be difficult to distinguish genuine from spurious matches for fragmented genes. Fragmentation is common among *K. pneumoniae* K loci, particularly when the genomes were sequenced using the Illumina technology with the Illumina XT library prep (see :ref:`FAQs<fragmented-Klebs-faq>` for more details).  



Untypeable loci
-----------------------

These are loci that do not meet the above criteria. **We recommend that users do not accept these results** unless
they are able to perform manual exploration of the data.

.. Problems:

Problems
=========
* ``?`` = the match was in a multiple pieces, possibly due to a poor match or discontiguous assembly. The number of pieces is indicated by the integer directly following the ``?`` symbol).
* ``-`` = genes expected in the locus were not found.
* ``+`` = extra genes were found in the locus.
* ``*`` = one or more expected genes was found but with translated identity below the minimum threshold.
* ``!`` = one or more genes was found but truncated


Exploring the other columns
=============================

Many users will not want or need to look beyond the columns described above. However, the rest of the Kaptive output may be useful for those wishing to investigate loci marked with ``Problems`` or explore locus variation in more detail. Interesting features include:

* Missing genes may indicate a novel locus or deletion variant is present in the assembly
* Length discrepencies can indicate a novel locus or deletion variant is present in the assembly. For *Klebsiella* K loci, positive length discrepencies of >700bp often indicate insertion sequence insertions resulting in so called 'IS variants' of the locus.   
* Other genes inside the locus may indicate a novel locus (with some exceptions, see the :ref:`FAQs<extra_genes_faq>`)
* Truncated genes may have an impact on the resultant phenotype. Kaptive will consider truncations when reporting predicting phenotypes, but it currenetly considers only gene truncations for which there is good supporting evidence in the literature, and such evidence is very limited.  


See the `tutorials <https://klebnet.org/training/>`_ for our tips on investigating loci in more detail outside of Kaptive.
