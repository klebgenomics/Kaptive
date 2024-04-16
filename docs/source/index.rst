.. Kaptive documentation master file, created by
   sphinx-quickstart on Tue Feb 27 16:31:58 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:

   Installation
   Usage
   Method
   Outputs
   Interpreting-the-results
   Databases
   FAQs

########################
Introducing Kaptive 3
########################

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.10252121.svg
    :target: https://doi.org/10.5281/zenodo.10252121
    :alt: DOI

.. image:: https://img.shields.io/github/stars/klebgenomics/Kaptive
   :alt: GitHub Repo stars

.. image:: https://readthedocs.org/projects/kaptive/badge/?version=latest
    :target: https://kaptive.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: License

.. image:: https://img.shields.io/conda/vn/bioconda/kaptive
   :target: https://anaconda.org/bioconda/kaptive
   :alt: Conda

About
==========
Kaptive is a system for surface polysaccharide typing from bacterial genome sequences. It consists of two main components:

#. Curated reference databases of surface polysaccharide gene clusters (loci).
#. A command-line interface (CLI) with three modes:

   -  **assembly**: surface polysaccharide typing from assemblies
   -  **extract**: extract features from Kaptive databases in different formats
   -  **convert**: convert Kaptive results to different formats

Kaptive can be found:

* On `GitHub <https://github.com/klebgenomics/Kaptive>`_ alongside the reference databases.
* On `PathogenWatch <https://pathogen.watch/>`_ where it is used to serotype *Klebsiella* and *Acinetobacter baumannii* isolates.
* As a GUI in `Kaptive Web <https://kaptive-web.erc.monash.edu/>`_ (`source code <https://github.com/kelwyres/Kaptive-Web>`_),
  which includes a `third-party database <https://github.com/aldertzomer/vibrio_parahaemolyticus_genomoserotyping>`_ for
  *Vibrio parahaemolyticus* [#vpara]_.

Citation
==========
* If you use Kaptive and the *Klebsiella* K or O locus databases in your research, please cite [#kaptive]_ and [#kaptive2]_.
* If you use `Kaptive Web <https://kaptive-web.erc.monash.edu/>`_, please cite [#kaptiveweb]_.
* If you use the *A. baumannii* K or OC locus database(s) in your research please cite [#aci]_ and [#aciupdate]_.
* Lists of papers describing each of the individual *A. baumannii* reference loci can be found `here <https://github.com/katholt/Kaptive/tree/master/extras>`_.

Tutorial
==========
Step-by-step `video <https://klebnet.org/training/>`_ and
`documented <https://docs.google.com/document/d/1aggwBCGu1CfsduOoKI0e6TRYOGtwwSceSBdKKkjCisA/edit?usp=sharing>`_
tutorials are available, covering:

* Kaptive's features and their scientific rationale
* How to run Kaptive
* Examples, illustrating how to run and interpret results
* Further investigations (e.g. exploring novel loci, IS insertions)

.. note::
 The tutorials are based on Kaptive 2.0, but the principles are similar for Kaptive 3.0.

People
==========
- `Dr. Tom Stanton <https://wyreslab.com/>`_
- `Dr. Kelly L. Wyres <https://wyreslab.com/>`_
- `Prof. Kathryn E. Holt <holtlab.net>`_
- `Dr. Ryan R. Wick <https://rrwick.github.io/>`_
- `A/Prof. Johanna Kenyon <https://experts.griffith.edu.au/45350-johanna-kenyon>`_

`Contact Kelly and Tom <kaptive.typing@gmail.com>`_ for help with Kaptive, or to report bugs or request features.


References
=============
.. [#aciupdate] Cahill, S.M., Hall, R.M., Kenyon, J.J., 2022. An update to the database for *Acinetobacter baumannii* capsular polysaccharide locus typing extends the extensive and diverse repertoire of genes found at and outside the K locus. *Microbial Genomics* 8. https://doi.org/10.1099/mgen.0.000878
.. [#kaptive2] Lam, M.M.C., Wick, R.R., Judd, L.M., Holt, K.E., Wyres, K.L., 2022. Kaptive 2.0: updated capsule and lipopolysaccharide locus typing for the *Klebsiella pneumoniae* species complex. *Microbial Genomics* 8. https://doi.org/10.1099/mgen.0.000800
.. [#vpara] Van Der Graaf-van Bloois, L., Chen, H., Wagenaar, J.A., Zomer, A.L., 2023. Development of Kaptive databases for *Vibrio parahaemolyticus* O- and K-antigen genotyping. *Microbial Genomics* 9. https://doi.org/10.1099/mgen.0.001007
.. [#kaptiveweb] Wick, R.R., Heinz, E., Holt, K.E., Wyres, K.L., 2018. Kaptive Web: User-Friendly Capsule and Lipopolysaccharide Serotype Prediction for *Klebsiella* Genomes. *Journal of Clinical Microbiology* 56, 10.1128/jcm.00197-18. https://doi.org/10.1128/jcm.00197-18
.. [#aci] Wyres, K.L., Cahill, S.M., Holt, K.E., Hall, R.M., Kenyon, J.J., 2020. Identification of *Acinetobacter baumannii* loci for capsular polysaccharide (KL) and lipooligosaccharide outer core (OCL) synthesis in genome assemblies using curated reference databases compatible with Kaptive. *Microbial Genomics* 6. https://doi.org/10.1099/mgen.0.000339
.. [#kaptive] Wyres, K.L., Wick, R.R., Gorrie, C., Jenney, A., Follador, R., Thomson, N.R., Holt, K.E., 2016. Identification of *Klebsiella* capsule synthesis loci from whole genome data. *Microbial Genomics* 2, e000102. https://doi.org/10.1099/mgen.0.000102
