====================================
Interpreting the results
====================================

Kaptive will indicate the best matching locus and types, and its confidence in the match:

* **Perfect** = the locus was found in a single piece with 100% coverage and 100% identity to the reference.

* **Very high** = the locus was found in a single piece with ≥99% coverage and ≥95% identity, with no missing genes and no extra genes.

* **High** = the locus was found in a single piece with ≥99% coverage, with ≤ 3 missing genes and no extra genes.

* **Good** = the locus was found in a single piece or with ≥95% coverage, with ≤ 3 missing genes and ≤ 1 extra genes.

* **Low** = the locus was found in a single piece or with ≥90% coverage, with ≤ 3 missing genes and ≤ 2 extra genes.

* **None** = did not qualify for any of the above.

We encourage users to make their own decisions about the confidence levels they are happy to accept for their specific
analyses, but as a general rule **we recommend users exclude matches with 'Low' or 'None' confidence** unless additional
investigations are made. Inspection of the additional columns in the output table can give some indications of the
causes of low confidence matches, but further exploration using BLASTn comparisons and visualisations are often required
e.g. using `Bandage assembly graph viewer <https://rrwick.github.io/Bandage/>`_ or the `Artemis Comparison Tool
<https://www.sanger.ac.uk/tool/artemis-comparison-tool-act/>`_. You can find examples in the `tutorials <https://klebnet.org/training/>`_.

**Assembly fragmentation is a common cause of low confidence matches among Klebsiella genomes.**
This can be identified by a `?` in the 'Problems' column, combined with low BLASTn coverage and usually one or more
missing genes. Unfortunately, the *Klebsiella* K-locus is a particularly difficult region of the genome to sequence due
to its low GC content compared to the rest of the chromosome, and this is particularly evident for draft genomes
assembled from Illumina data generated with the Nextera XT library prep kit. These read sets may have little or no
read coverage of parts of the K-locus, meaning it is impossible to fully assemble it. Since Kaptive uses BLASTn coverage
as a key metric in its confidence scoring, it will generally assign low confidence scores to fragmented loci,
particularly if parts of the locus are totally absent from the assembly. In these cases, Kaptive may also incorrectly
preference the shortest locus in the database as the best match - KL107 - because it has the highest BLASTn coverage.
**Therefore, 'Low' and 'None' confidence matches to KL107 should always be excluded**. However, 'Low' and 'None'
confidence hits to other loci may be more reliable. See the `FAQs <https://github.com/klebgenomics/Kaptive/wiki/FAQs#how-should-i-report-low-or-none-confidence-matches-for-klebsiella-genomes>`_ for our recommendations on dealing with these.

Below are some examples of Kaptive results and suggested interpretations. The gene details columns of the table have
been excluded for brevity, as they can be quite long.

'Perfect' and 'Very high' confidence matches
==============================================

Example 1
-----------

============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
Assembly     Best match locus   Best match type   Confidence   Problems   Coverage   Identity   Length discrepancy   Expected genes in locus   Expected genes in locus, details   Missing expected genes   Other genes in locus   Other genes in locus, details   Expected genes outside locus   Expected genes outside locus, details   Other genes outside locus   Other genes outside locus, details
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
assembly_1   KL1                K1                Perfect                 100%       100%       0 bp                 20 / 20 (100%)            ...                                                         0                                                      0                                                                      2                           ...
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================

This is a case where our assembly perfectly matches a known K locus, KL1. There are no characters in the 'Problems' column, the coverage and identity are both 100%, with no length discrepancy. A couple of other low-identity K locus genes hits were elsewhere in the assembly, but that's not abnormal and no cause for concern. The KL1 locus is known to encode the K1 serotype, as indicated in the 'Best match type' column.

Example 2
-----------

============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
Assembly     Best match locus   Best match type   Confidence   Problems   Coverage   Identity   Length discrepancy   Expected genes in locus   Expected genes in locus, details   Missing expected genes   Other genes in locus   Other genes in locus, details   Expected genes outside locus   Expected genes outside locus, details   Other genes outside locus   Other genes outside locus, details
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
assembly_2   KL1                K1                Very high               99.94%     99.81%     -22 bp               20 / 20 (100%)            ...                                                         0                                                      0                                                                      2                           ...
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================

This is a case where our assembly very closely matches a known K locus, KL1. There are no characters in the 'Problems' column, the coverage and identity are both high (although not 100%), the length discrepancy is low, and all expected genes were found with high identity. A couple of other low-identity K locus genes hits were elsewhere in the assembly, but that's not abnormal and no cause for concern.

Overall, this is a nice, solid match for the KL1 locus, which is known to encode the K1 serotype, as indicated in the 'Best match type' column.

Example 3
-----------

============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
Assembly     Best match locus   Best match type   Confidence   Problems   Coverage   Identity   Length discrepancy   Expected genes in locus   Expected genes in locus, details   Missing expected genes   Other genes in locus   Other genes in locus, details   Expected genes outside locus   Expected genes outside locus, details   Other genes outside locus   Other genes outside locus, details
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
assembly_3   KL102              unknown           Very high    *          99.84%     95.32%     +97 bp               19 / 19 (100%)            ...                                                         0                                                      0                                                                      2                           ...
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================

This case shows an assembly that matches the KL102 locus sequence (encoding an unknown serotype), but not as closely as our previous examples. The `*` character indicates that one or more of the expected genes falls below the identity threshold (default 95%). The 'Expected genes in K locus, details' columns, excluded here for brevity, would show the identity for each gene.

Our sample still almost certainly has a K locus of KL102, but it has diverged a bit more from the reference, possibly due to mutation and/or recombination. 

'High' and 'Good' confidence matches
======================================

Example 4
-----------

============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
Assembly     Best match locus   Best match type   Confidence   Problems   Coverage   Identity   Length discrepancy   Expected genes in locus   Expected genes in locus, details   Missing expected genes   Other genes in locus   Other genes in locus, details   Expected genes outside locus   Expected genes outside locus, details   Other genes outside locus   Other genes outside locus, details
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
assembly_4   KL38               K38               High         -          99.99%     97.45%     -1bp                 17 / 18 (94.4%)           ...                                KL38-CDS12-rmlD          0                                                      0                                                                      1                           ...
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================

Here the assembly has both a high coverage and high identity match to the KL38 locus, which is known to encode the K38 serotype (as indicated in the 'Best match type' column). The locus is found in a single piece in the assembly but with a length discrepancy of -1bp. The KL38-CDS12-rmlD gene is also missing. This may seem confusing given the very high coverage BLAST match to the reference locus, but it is worth remembering that the locus coverage is calculated using BLASTn, whereas the genes are detected using tBLASTn. So in this case, we may suspect that there is a frameshift mutation within KL38-CDS12-rmlD, resulting in a predicted protein truncation that drops the gene coverage below the threshold for detection (default 90%). Such a mutation might have an important impact on serotype, and so you may choose to investigate further.

Example 5
-----------

============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
Assembly     Best match locus   Best match type   Confidence   Problems   Coverage   Identity   Length discrepancy   Expected genes in locus   Expected genes in locus, details   Missing expected genes   Other genes in locus   Other genes in locus, details   Expected genes outside locus   Expected genes outside locus, details   Other genes outside locus   Other genes outside locus, details
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
assembly_5   KL2                K2                Good         ?-         99.95%     98.38%     n/a                  17 / 18 (94.4%)           ...                                K2-CDS17-manB            0                                                      0                                                                      1                           ...
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================

Here is a case where our assembly matched a known K locus well (high coverage and identity) but with a couple of problems. First, the `?` character indicates that the K locus sequence was not found in one piece in the assembly. Second, one of the expected genes (K2-CDS17-manB) was not found in the gene BLAST search.

In cases like this, it may be worth examining in more detail outside of Kaptive. For this example, such an examination may have revealed that the assembly was poor (broken into many small pieces) and the *manB* gene happened to be split between two contigs. So the *manB* gene isn't really missing, it's just broken in two. If that were the case, the sample most likely is a very good match for KL2, but the poor assembly quality made it difficult for Kaptive to determine that automatically.

Alternatively, this example may represent an insertion sequence variant of the KL2 locus e.g. if an insertion sequence (a mobile genetic element) has inserted within *manB* causing a gene interruption. Since insertion sequences are often present in multiple copies in a genome, they usually cause assembly fragmentation. 

Depending on your specific analysis question, you may be happy to simply assign this assembly as a match to KL2 (e.g. if your question is around K locus epidemiology), or you may wish to distinguish it is an insertion sequence variant (e.g. if your goal is to predict serotypes, for which the impact of such an insertion may or may not be known). 

'Low' and 'None' confidence matches
=====================================

Example 6
-----------

============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
Assembly     Best match locus   Best match type   Confidence   Problems   Coverage   Identity   Length discrepancy   Expected genes in locus   Expected genes in locus, details   Missing expected genes   Other genes in locus   Other genes in locus, details   Expected genes outside locus   Expected genes outside locus, details   Other genes outside locus   Other genes outside locus, details
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
assembly_6   KL116              unknown           Low          ?-         92.26%     99.98%     n/a                  18 / 21 (85.7%)           ...                                ...                      0                                                      0                                                                      2                           ...
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================

In this example the assembly has a reasonable coverage match to the KL116 locus with high identity, but 3 genes are missing and the locus was not found in a single piece in the assembly (`?` in the 'Problems' column). This may be because the assembly is fragmented due to sequencing depth and/or assembly problems, or it may be because the assembly doesn't carry a true match to the KL116 locus. There is no way to distinguish these possibilities without further investigation outside of Kaptive. 


Example 7
-----------

============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
Assembly     Best match locus   Best match type   Confidence   Problems   Coverage   Identity   Length discrepancy   Expected genes in locus   Expected genes in locus, details   Missing expected genes   Other genes in locus   Other genes in locus, details   Expected genes outside locus   Expected genes outside locus, details   Other genes outside locus   Other genes outside locus, details
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================
assembly_7   KL137              unknown           None         ?-*        77.94%     83.60%     n/a                  15 / 20 (75%)             ...                                ...                      0                                                      0                                                                      5                           ...
============ ================== ================= ============ ========== ========== ========== ==================== ========================= ================================== ======================== ====================== =============================== ============================== ======================================= =========================== ====================================

In this case, Kaptive did not find a close match to any known K locus sequence. The best match was to KL137, but BLAST only found alignments for 78% of the KL137 sequence, and only at 84% nucleotide identity. Five of the twenty KL137 genes were not found, and the 15 which were found had low identity. The assembly sequences matching KL137 did not come in one piece (indicated by `?`), possibly due to assembly problems, but more likely due to the fact that our sample is not in fact KL137 but rather has some novel K locus that was not in our reference inputs.

A case such as this should not be reported as KL137 unless closer examination is completed outside of Kaptive. It could be a novel locus, and you may wish to extract and annotate the locus sequence from the assembly. 

**If you do believe that you've found a novel locus, please consider submitting for inclusion in the relevant database.**
