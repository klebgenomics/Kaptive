====================================
How does Kaptive work?
====================================

Given a novel genome and a database of known loci (K, O or OC), Kaptive will help a user to decide whether their sample has a known or novel locus. It carries out the following for each input assembly:

* BLAST for all known locus nucleotide sequences (using `blastn`) to identify the best match ('best' defined as having the highest coverage).

* Extract the region(s) of the assembly which correspond to the BLAST hits (i.e. the locus sequence in the assembly) and save it to a FASTA file.

* BLAST for all known locus genes (using `tblastn`) to identify which expected genes (genes in the best matching locus) are present/missing and whether any unexpected genes (genes from other loci) are present.

* Output a summary to a table file.

In cases where your input assembly closely matches a known locus, Kaptive should make that obvious. When your assembly has a novel type, that too should be clear. However, Kaptive cannot reliably extract or annotate locus sequences for totally novel types – if it indicates a novel locus is present then extracting and annotating the sequence is up to you! Very poor assemblies can confound the results, so be sure to closely examine any case where the locus sequence in your assembly is broken into multiple pieces.

If you think you have found a novel locus that should be added to one of the databases distributed with Kaptive please `contact us <mailto:kaptive.typing@gmail.com>`_.
