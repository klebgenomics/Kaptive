
**************************************
Interpreting the results
**************************************

The four most important columns in the Kaptive :ref:`tabular <Tabular>` output:

* ``Best match locus`` indicating the locus that is best supported by teh avaialble sequence data
* ``Best match type`` indicating the predicted phenotype based on that associated with the `Best match locus` and taking into account any special phenotype logic e.g. where a known essential genes is truncated, Kaptive will report the `Best match type` as 'Capsule null'. 
* ``Match confidence`` indicates a qualitative level of confidence that the reported `Best match locus` is correct (see below).
* ``Problems`` indicates the features of the assembly locus that may impact ``Match confidence``

.. note::
  The ``Match confidence`` relates to the locus match, and does not a direct indication of confidence in the ``Best match type``.  



Confidence score  
=================
Kaptive will indicate the best matching locus and its confidence in the locus match.


Good confidence loci  
---------------------

The locus was found in a single piece in the query assembly with:

* no missing genes (as determined by at least 50% nucletide mapping coverage)  
* no unexpected genes (genes from other loci) inside the locus region of the assembly  
* no genes below the minimum translated identity threshold (82.5% for *K. pneumoniae*, 85% for *A. baumanii*)  
* at least 90% overall nucletide identity.  

OR

The locus was found in more than one piece in the query assembly with:

* no more than 1 missing gene  
* no more than 1 unespected gene (genes from other loci) inside the locus region of the assembly  
* no genes below the minimum translated identity threshold  
* at least 90% overall nucletide identity.

These criteria were designed in consideration of the locus definition rules (i.e. that each locus represents a unique set of genes defined at a givin minimum translated identity threshold) and following systematic analysis of Kaptive outputs for draft genome assemblies compared against manually confirmed loci determined from matched completed genomes.

We allow some flexibility with regards to missing genes or additional genes found within the locus when this region of the query assembly is fragmented, because it can be difficult to distinguish genuine from spurious matches for fragmented genes. Fragmentation is common among *K. pneumoanie* K-loci, particularly when the genomes were sequenced using the Illumina technology with the Illumina XT library prep (see :ref:`FAQs<fragmented-Klebs-faq>` for more details).  


Low confidence loci
-----------------------

The locus was found in a single piece in the query assembly with: 

* no missing genes (as determined by at least 50% nucletide mapping coverage)  
* no unexpected genes (genes from other loci) inside the locus region of the assembly  
* no genes below the minimum translated identity threshold (82.5% for *K. pneumoniae*, 85% for *A. baumanii*)  
* less than 90% overall nucletide identity.  

OR

The locus was found in more than one piece in the query assembly with:

* no more than 1 missing gene  
* no more than 1 unespected gene (genes from other loci) inside the locus region of the assembly  
* no genes below the minimum translated identity threshold  
* less than 90% overall nucletide identity.  

These criteria are essentially the same as those for ``good`` matches, with exception that overall nucleotide identity is below 90%. This reflects our observations from systematic testing where correctly identified loci generally share >90% identity with the reference locus, while novel loci generally had lower identity. However, based on the locus definition rules (where minimum gene identities are set at 82-85%), it is technically possible for an assembly to carry a locus matching the reference with <90% nucleotide identity. 

**We therefore recommend that users investigate low confidence hits** through inspection of the rest of the Kaptive output and/or investiogations with other tools. For well characterised loci such as the *K. pneumoniae* K/O loci and *A. baumanii* K/OC loci, we expect only a minority of genomes to return a ``low`` confidence match.


Untypeable loci
-----------------------

These are loci that do not meet the above criteria. **We recommend that users do not accept these results** unless they are able to perform manual exploration of the data.


Problems
=========
* ``?`` = the match was in a multiple pieces, possibly due to a poor match or discontiguous assembly. The number of pieces is indicated by the integer directly following the `?` symbol).
* ``-`` = genes expected in the locus were not found.
* ``+`` = extra genes were found in the locus.
* ``*`` = one or more expected genes was found but with translated identity below the minimum threshold.
* ``!`` = one or more genes was found but truncated


Exploring the other columns
=============================

Moany users will not want or need to look beyond the columns described above. However, the rest of the Kaptive output may be useful for those wishing to investigate ``low`` confidence matches or explore locus variation in more detail. Interesting features include:

* Missing genes may indicate a novel locus or deletion variant is present in the assembly
* Length discrepencies can indicate a novel locus or deletion variant is present in the assembly. For *Klebsiella* K loci, positive length discrepencies of >700bp often indicate insertion sequence insertions resulting in so called 'IS variants' of the locus.   
* Other genes inside the locus may indicate a novel locus (with some exceptions, see the :ref:`FAQs<extra_genes_faq>`)
* Truncated genes may have an impact on the resultant phenotype. Kaptive will condier truncations when reporting predicting phenotypes, but it currenetly considers only gene truncations for which there is good supprting evidence in the literature, and such evidence is very limited.  


See the `tutorials <https://klebnet.org/training/>`_ for our tips on investigating loci in more detail outside of Kaptive.

