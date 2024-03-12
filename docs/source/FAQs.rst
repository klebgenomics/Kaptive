====================================
FAQs
====================================

Why do I have a large number of *Klebsiella* genomes with 'Low' or 'None' confidence matches
================================================================================================

Unfortunately, the *Klebsiella* K-locus is a particularly difficult region of the genome to sequence due to its low GC content compared to the rest of the chromosome, which can result in assembly fragmentation. This is particularly evident for draft genomes assembled from Illumina data generated with the Nextera XT library prep kit. These read sets may have little or no read coverage of parts of the K-locus, meaning it is impossible to fully assemble it. Since Kaptive uses BLASTn coverage as a key metric in its confidence scoring, it will generally assign low confidence scores to fragmented loci, particularly if parts of the locus are totally absent from the assembly. In these cases, Kaptive may also incorrectly preference the shortest locus in the database as the best match - KL107 - because it has the highest BLASTn coverage. These 'Low' or 'None' confidence matches to KL107 should be excluded from analyses and treated as 'untypeable'. 

How should I report 'Low' or 'None' confidence matches for *Klebsiella* genomes
===================================================================================

The most conservative and simplest approach is to treat all of these genomes as 'untypeable'. Alternatively, you could undertake some manual investigations to figure out what is causing the 'Low' and 'None' confidence calls (see the `tutorials <https://klebnet.org/training/>`_ for our tips). However, manual investigations may not be possible for genome collections with a high number of 'Low' and 'None' confidence calls (e.g. due to assembly fragmentation, explained above). As an alternative, you may choose to accept some of the 'Low' or 'None' confidence calls without further investigation if: 

* there is `?` in the 'Problems' column (indicating that the K-locus is fragmented in the assembly) and;

* there are no additional genes reported within the locus (excluding *ugd* - see below) and;  

* the Best Match is not reported as KL107 - which is often incorrectly identified among highly fragmented assemblies because it is the shortest locus in the reference database (and therefore returns a high BLASTn coverage).


Why are there locus genes found outside the locus?
======================================================

For *Klebsiella* K loci in particular, a number of the K-locus genes are orthologous to genes outside of the K-locus region of the genome. E.g the *Klebsiella* K-locus *man* and *rml* genes have orthologues in the LPS (lipopolysacharide) locus; so it is not unusual to find a small number of genes "outside" the locus.

However, if you have a large number of genes (>5) outside the locus it may mean that there is a problem with the locus match, or that your assembly is very fragmented or contaminated (contains more than one sample).

How can my sample be missing locus genes when it has a full-length, high identity locus match?
==================================================================================================

Kaptive uses 'tblastn' to screen for the presence of each locus gene with a coverage threshold of 90% (default). A single non-sense mutation or small indel in the centre of a gene will interrupt the 'tblastn' match and cause it to fall below the 90% threshold. However, such a small change has only a minor effect on the nucleotide 'blast' match across the full locus.

Why does the *Klebsiella* K-locus region of my sample contain a *ugd* gene matching another locus?
===========================================================================================================

A small number of the original *Klebsiella* K locus references are truncated, containing only a partial *ugd* sequence. The reference annotations for these loci do not include *ugd*, so are not identified by the 'tblastn' search. Instead <b>Kaptive</b> reports the closest match to the partial sequence (if it exceeds the 90% coverage threshold). 

Why has the best matching locus changed after I reran my analysis with an updated version of the database?
===================================================================================================================

The databases are updated as novel loci are discovered and curated. If your previous match had a confidence call of 'Low' or 'None' but your new match has higher confidence, this indicates that your genome contains a locus that was absent in the older version of the database. So nothing to worry about here.

But what if your old match and your new match have 'Good' or better confidence levels?

If your old match had 'Perfect' or 'Very High' confidence, please post an issue to the issues page, as this may indicate a problem with the new database!

If your old match had 'Good' or 'High' confidence please read on...

Polysaccharide loci are subject to frequent recombinations and rearrangements, which generates new variants. As a result, a small number of pairs of loci share large regions of homology e.g. the *Klebsiella* K-locus KL170 is very similar to KL101, and in fact seems to be a hybrid of KL101 plus a small region from KL106. 

Kaptive can accurately distinguish the KL101 and KL170 loci when it is working with high quality genome assemblies, but this task is much trickier if the assembly is fragmented. This means that matches to KL101 that were reported using an early version of the K-locus database might be reported as KL170 when using a later version of the database. However, this should only occur in instances where the K-locus is fragmented in the genome assembly and in that case Kaptive will have indicated 'problems' with the matches (e.g. '?' indicating fragmented assembly or '-' indicating that an expected gene is missing), and the corresponding confidence level will be at the lower end of the scale (i.e. 'Good' or 'High', but not 'Very High' or 'Perfect').

You may want to confirm the correct locus manually, e.g. using `Bandage <https://rrwick.github.io/Bandage/>`_ to BLAST the corresponding loci in your genome assembly graph. 
