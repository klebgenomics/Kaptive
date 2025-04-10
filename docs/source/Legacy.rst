*******************************************************
Legacy Database Information
*******************************************************

.. note::
    This page contains information about legacy Kaptive databases for backwards compatibility and interpretation of
    old results.

.. _Legacy-Klebsiella-O-locus-database:

Legacy *Klebsiella* O locus database
------------------------------------------

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

========= ==================== ================================== ================================= ======================================== ==================================
O locus   Extra genes          Kaptive < v2.0 (locus\ :sup:`a`)   Kaptive v2.0+ (locus\ :sup:`a`)   Kaptive v2.0 - v2.0.7 (type\ :sup:`b`)   Kaptive v2.0.8+ (type\ :sup:`b`)
========= ==================== ================================== ================================= ======================================== ==================================
O1/O2v1   none                 O2v1                               O1/O2v1                           O2a                                      O2a
O1/O2v2   none                 O2v2                               O1/O2v2                           O2afg                                    O2afg
O1/O2v3   none                 Na                                 O1/O2v3                           O2a                                      O2a
O1/O2v1   *wbbYZ*              O1v1                               O1/O2v1                           Na                                       O1ab
O1/O2v2   *wbbYZ*              O1v2                               O1/O2v2                           Na                                       O1ab
O1/O2v3   *wbbYZ*.             Na                                 O1/O2v3                           Na                                       O1ab
O1/O2v1   *wbbY* only          O1v1                               O1/O2v1                           O1                                       O1a
O1/O2v2   *wbbY* only          O1v2                               O1/O2v2                           O1                                       O1a
O1/O2v3   *wbbY* only          Na                                 O1/O2v3                           O1                                       O1a
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
========= ==================== ================================== ================================= ======================================== ==================================


- :sup:`a` as reported in the ‘Best match locus’ column in the Kaptive output.
- :sup:`b` predicted antigenic serotype reported in the 'Best match type' column in the Kaptive output (v2.0 and above).
- Na- not applicable