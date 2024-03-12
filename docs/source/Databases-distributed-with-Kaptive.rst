====================================
Databases distributed with Kaptive
====================================

Kaptive is distributed with databases for detection of *Klebsiella pneumoniae* species complex and *Acinetobacter baumanii* surface antigen synthesis loci in the `reference_database <https://github.com/katholt/Kaptive/tree/master/reference_database>`_ directory, (see details below). You can also generate your own databases for use with Kaptive by following these guidelines.

The existing databases were developed and curated by `Kelly Wyres <https://holtlab.net/kelly-wyres/>`_ (*Klebsiella*)
and `Johanna Kenyon <https://research.qut.edu.au/infectionandimmunity/projects/bacterial-polysaccharide-research/>`_
(*A. baumannii*).

A third-party Kaptive database is available for *Vibrio parahaemolyticus* `K and O loci <https://github.com/aldertzomer/vibrio_parahaemolyticus_genomoserotyping>`_,
created by Aldert Zomer and team (see `preprint <https://doi.org/10.1101/2021.07.06.451262>`_).
The database can be `downloaded <https://github.com/aldertzomer/vibrio_parahaemolyticus_genomoserotyping>`_ and
used as input to command-line Kaptive, it is also available in the online tool `Kaptive-Web <https://kaptive-web.erc.monash.edu/>`_
along with our *Klebsiella* and *A. baumannii* databases.

We are always keen to expand the utility of Kaptive for the research community, so if you have created a database that \
you feel will be useful for others and you are willing to share this resource, please get in touch via the
`issues page <https://github.com/katholt/Kaptive/issues>`_ or `email <mailto:kaptive.typing@gmail.com>`_.

Similarly, if you have identified new locus variants not currently in the existing databases, please let us know!


*Klebsiella* K locus databases
================================

The *Klebsiella* K locus primary reference database (``Klebsiella_k_locus_primary_reference.gbk``) comprises full-length
(*galF* to *ugd*) annotated sequences for each distinct *Klebsiella* K locus, where available:

* KL1 - KL77 correspond to the loci associated with each of the 77 serologically defined K-type references, for which
  the corresponding predicted serotypes are K1-K77, respectively.
* KL101 and above are defined from DNA sequence data on the basis of gene content, and are not currently associated with
  any defined serotypes.

Note that insertion sequences (IS) are excluded from this database since we assume that the ancestral sequence was
likely IS-free and IS transposase genes are not specific to the K locus.

Synthetic IS-free K locus sequences were generated for K loci for which no naturally occurring IS-free variants have
been identified to date.

The variants database (``Klebsiella_k_locus_variant_reference.gbk``) comprises full-length annotated sequences for
variants of the distinct loci:

* IS variants are named as KLN -1, -2 etc e.g. KL15-1 is an IS variant of KL15.
* Deletion variants are named KLN-D1, -D2 etc e.g. KL15-D1 is a deletion variant of KL15.

Note that KL156-D1 is included in the primary reference database since no full-length version of this locus has been
identified to date.

We recommend screening your data with the primary reference database first to find the best-matching K locus. If you
have poor matches or are particularly interested in detecting variant loci you should try the variant database.

**WARNING**: If you use the variant database please inspect your results carefully and decide for yourself what
constitutes a confident match! Kaptive is not optimised for accurate variant detection.

Database versions:

* Kaptive releases v0.5.1 and below include the original *Klebsiella* K locus databases, as described in
  `Wyres, K. et al. Microbial Genomics 2016. <http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000102>`_
* Kaptive v0.6.0 and above include four novel primary *Klebsiella* K locus references defined on the basis of gene
  content (KL162-KL165) in `Wyres et al. Genome Medicine 2020 <https://pubmed.ncbi.nlm.nih.gov/31948471/>`_.
* Kaptive v0.7.1 and above contain updated versions of the KL53 and KL126 loci (see table below for details).
  The updated KL126 locus sequence is described in `McDougall, F. et al. Research in Microbiology 2021 <https://pubmed.ncbi.nlm.nih.gov/34506927/>`_.
* Kaptive v0.7.2 and above include a novel primary *Klebsiella* K locus reference defined on the basis of gene content
  (KL166), described in `Le, MN. et al. Microbial Genomics 2022 <https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000827>`_.
* Kaptive v0.7.3 and above include four novel primary *Klebsiella* K locus references defined on the basis of gene
  content (KL167-KL170), described in `Gorrie, C. et al. Nature Communications 2022. <https://www.nature.com/articles/s41467-022-30717-6>`_
* Kaptive v2.0 and above include 16 novel primary *Klebsiella* K locus references defined on the basis of gene content
  (KL171-KL186) and described in `Lam, M.M.C et al. Microbial Genomics 2022. <https://doi.org/10.1099/mgen.0.000800>`_

Changes to the *Klebsiella* K locus primary reference database:

======= =================================================================================================== ========================================================================== ================ =====================
Locus   Change                                                                                              Reason                                                                     Date of change   Kaptive version no.
======= =================================================================================================== ========================================================================== ================ =====================
KL53    Annotation update: *wcaJ* changed to *wbaP*                                                         Error in original annotation                                               21 July 2020     v 0.7.1
KL126   Sequence update: new sequence from isolate FF923 includes *rmlBADC* genes between *gnd* and *ugd*   Assembly scaffolding error in original sequence from isolate A-003-I-a-1   21 July 2020     v 0.7.1
======= =================================================================================================== ========================================================================== ================ =====================


*Klebsiella* O locus database
===============================

The *Klebsiella* O locus database (``Klebsiella_o_locus_primary_reference.gbk``) contains annotated sequences for 13
distinct *Klebsiella* O loci.

O locus classification requires some special logic, as the O1 and O2 serotypes are associated with the same loci and
the distinction between O1 and each of the four defined O2 subtypes (O2a, O2afg, O2ac, O2aeh) is determined by the
presence/absence of 'extra genes' elsewhere in the chromosome as indicated in the table below. Kaptive therefore looks
for these genes to predict antigen (sub)types. (Note that the original implementation of O locus typing in Kaptive
(< v2.0) distinguished O1 and O2 but not the O2 subtypes.)

Read more about the O locus and its classification here: `The diversity of *Klebsiella* pneumoniae surface polysaccharides <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320592/>`_.

Find out about the genetic determinants of O1 and O2 (sub)types here: `Molecular basis for the structural diversity in serogroup O2-antigen polysaccharides in *Klebsiella pneumoniae* <https://pubmed.ncbi.nlm.nih.gov/29602878/>`_.

Find out about the O1 glycoforms and their genetic determinants here: `Identification of a second glycoform of the clinically prevalent O1 antigen from *Klebsiella pneumoniae*. <https://doi.org/10.1073/pnas.2301302120>`_

Database versions:

* Kaptive v0.4.0 and above include the original version of the *Klebsiella* O locus database, as described in `Wick, R. et al. J Clin Microbiol 2019 <http://jcm.asm.org/content/56/6/e00197-18>`_.
* Kaptive v2.0 and above include a novel O locus reference (O1/O2v3) and updated 'Extra genes' for prediction of O1 and O2 antigen (sub)types, as shown in the table below and described in `Lam, M.M.C et al. 2021. Microbial Genomics 2022. <https://doi.org/10.1099/mgen.0.000800>`_
* Kaptive v2.0.8 and above include:

  i) updated 'Extra genes' logic for prediction of O1 glycoforms, reported as O1a (isolate predicted to produce O1a
     only) and O1ab (isolate predicted to be able to produce both O1a and O1b glycoforms);

  ii) OL101 re-assigned as OL13 and its associated phenotype prediction updated to O13, to reflect the
      `description of the novel O13 polysaccharide structure <https://www.sciencedirect.com/science/article/pii/S0144861723010469>`_.

Genetic determinants of O1 and O2 outer LPS antigens as reported in Kaptive:

========= ==================== ================================== ================================= ======================================== ================================== ==
O locus   Extra genes          Kaptive < v2.0 (locus\ :sup:`a`)   Kaptive v2.0+ (locus\ :sup:`a`)   Kaptive v2.0 - v2.0.7 (type\ :sup:`b`)   Kaptive v2.0.8+ (type\ :sup:`b`)
========= ==================== ================================== ================================= ======================================== ================================== ==
O1/O2v1   none                 O2v1                               O1/O2v1                           O2a                                      O2a
O1/O2v2   none                 O2v2                               O1/O2v2                           O2afg                                    O2afg
O1/O2v3   none                 Na                                 O1/O2v3                           O2a                                      O2a
O1/O2v1   *wbbYZ*              O1v1                               O1/O2v1                           Na                                       O1ab
O1/O2v2   *wbbYZ*              O1v2                               O1/O2v2                           Na                                       O1ab
O1/O2v3   *wbbYZ*.             Na                                 O1/O2v3                           Na                                       O1ab
O1/O2v1   *wbbY* only          O1v1                               O1/O2v1                           O1                                       O1a
O1/O2v2   *wbbY* only          O1v2                               O1/O2v2                           O1                                       O1a
O1/O2v3   *wbbY* only          Na.                                O1/O2v3                           O1                                       O1a
O1/O2v1   *wbbY* OR *wbbZ*     O1/O2v1                            Na                                Na                                       Na
O1/O2v2   *wbbY* OR *wbbZ*     O1/O2v2                            Na                                Na                                       Na
O1/O2v3   *wbbY* OR *wbbZ*     Na                                 Na                                Na                                       Na
O1/O2v1   *wbmVW*              Na                                 O1/O2v1                           O2ac                                     O2ac
O1/O2v2   *wbmVW*              Na                                 O1/O2v2                           O2ac                                     O2ac
O1/O2v3   *wbmVW*              Na                                 O1/O2v3                           O2ac                                     O2ac
O1/O2v1   *gmlABD*             Na                                 O1/O2v1                           O2aeh                                    O2aeh
O1/O2v2   *gmlABD*             Na                                 O1/O2v2                           O2aeh                                    O2aeh
O1/O2v3   *gmlABD*             Na                                 O1/O2v3                           O2aeh                                    O2aeh
O1/O2v1   *wbbY* AND *wbmVW*   Na                                 O1/O2v1                           O1 (O2ac)\ :sup:`b`                      O1 (O2ac)\ :sup:`b`
O1/O2v2   *wbbY* AND *wbmVW*   Na                                 O1/O2v2                           O1 (O2ac)\ :sup:`b`                      O1 (O2ac)\ :sup:`b`
O1/O2v3   *wbbY* AND *wbmVW*   Na                                 O1/O2v3                           O1 (O2ac)\ :sup:`b`                      O1 (O2ac)\ :sup:`b`
========= ==================== ================================== ================================= ======================================== ================================== ==


- :sup:`a` as reported in the ‘Best match locus’ column in the Kaptive output.
- :sup:`b` predicted antigenic serotype reported in the 'Best match type' column in the Kaptive output (v2.0 and above).
- Na- not applicable

*Acinetobacter baunannii* K and OC locus databases
====================================================

The *A. baumannii* K (capsule) locus reference database (`Acinetobacter*baumannii*k*locus*primary_reference.gbk`)
contains annotated sequences for 241 distinct K loci.

The *A. baumannii* OC (lipooligosaccharide outer core) locus reference database (`Acinetobacter*baumannii*OC*locus*primary_reference.gbk`)
contains annotated sequences for 22 distinct OC loci.

**WARNING:** These databases have been developed and tested specifically for *A. baumannii* and may not be suitable for
screening other *Acinetobacter* species. You can check that your assembly is a true *A. baumannii* by screening for the
*oxaAB* gene e.g. using blastn.

Database versions:

* Kaptive v0.7.0 and above include the original *A. baumannii* K and OC locus databases, as described in
  `Wyres, KL. et al. Microbial Genomics 2020 <https://doi.org/10.1099/mgen.0.000339>`_.

* Kaptive v2.0.1 and above include 149 novel primary *A. baumannii* K locus references as described in
  Cahill, S.M. et al. 2022. An update to the database for *Acinetobacter baumannii* capsular polysaccharide locus typing
  extends the extensive and diverse repertoire of genes found at and outside the K locus.
  `Microbial Genomics <https://doi.org/10.1099/mgen.0.000878>`_.

* Kaptive v2.0.2 and above include special logic parameters that enable prediction of the capsule polysaccharide type
  based on KL or the detected combination of a specific KL with 'extra genes' elsewhere in the chromosome as indicated
  in the table below and described in Cahill, S.M. et al. 2022. An update to the database for *A. baumannii* capsular
  polysaccharide locus typing extends the extensive and diverse repertoire of genes found at and outside the K locus.
  `Microbial Genomics <https://doi.org/10.1099/mgen.0.000878>`_.

* Kaptive v2.0.5 and above includes a further 10 *A. baumannii* OC locus references (OCL13-OCL22) as described in
  Sorbello, B. et al. Identification of further variation at the lipooligosaccharide outer core locus in
  *Acinetobacter baumannii* genomes and extension of the OCL reference sequence database for Kaptive. *In prep*.


========= ================== ===================
K locus   Extra genes        K type
========= ================== ===================
KL1       None               K1
KL1       wzy-GI1,atr25-GI   K1-Wzy/Atr25-GI1
KL1       wzy-GI1            K1-Wzy-GI1
KL1       atr25-GI1          K1-Atr25-GI
KL5       None               K5
KL5       atr30-Ph           K5-Atr30-Ph
KL19      None               unknown (KL19)
KL19      wzy-GI1,atr25-GI   K19-Wzy/Atr25-GI1
KL19      wzy-GI1            K19-Wzy-GI1
KL24      None               unknown (KL24)
KL24      wzy-GI2            K24-Wzy-GI2
KL46      None               K46
KL46      atr29-Ph           K46-Atr29-Ph
KL127     None               K127
KL127     wzy-Ph1            K127-Wzy-Ph1
========= ================== ===================
