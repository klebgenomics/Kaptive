<img src="docs/assets/logo.png" alt="Kaptive" width="600"></p>

![GitHub Repo stars](https://img.shields.io/github/stars/klebgenomics/Kaptive)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/kaptive)
[![Conda](https://img.shields.io/conda/vn/bioconda/kaptive)](https://anaconda.org/bioconda/kaptive)

## About

Kaptive is a system for surface polysaccharide typing from bacterial
genome sequences. It consists of two main components:

1.  Curated reference [databases](docs/Databases.md#distributed-databases) of
    surface polysaccharide gene clusters (loci).
2.  A command-line interface (CLI) with three modes:
    - **assembly**: surface polysaccharide typing from assemblies
    - **extract**: extract features from Kaptive databases in different
      formats
    - **convert**: convert Kaptive results to different formats

Kaptive can be found:

- On [GitHub](https://github.com/klebgenomics/Kaptive) alongside the
  reference [databases](docs/Databases.md#distributed-databases).
- On [PathogenWatch](https://pathogen.watch/) where it is used to
  serotype *Klebsiella* and *Acinetobacter baumannii* isolates.
- As a GUI in [Kaptive Web](https://kaptive-web.erc.monash.edu/)
  ([source code](https://github.com/kelwyres/Kaptive-Web)), which
  includes a [third-party
  database](https://github.com/aldertzomer/vibrio_parahaemolyticus_genomoserotyping)
  for *Vibrio parahaemolyticus*[^1].

For full **documentation**, including install and usage instructions, **click [here](https://klebgenomics.github.io/Kaptive)**.

## Citation

- If you use Kaptive in your research, please cite[^2].
- If you use the *Klebsiella* K or O locus databases, please cite[^3]
  and[^4].
- If you use [Kaptive Web](https://kaptive-web.erc.monash.edu/), please
  cite[^5].
- If you use the *A. baumannii* K or OC locus database(s) in your
  research please cite[^6] and[^7].

## Tutorial

Step-by-step [video](https://klebnet.org/training/) and
[documented](https://docs.google.com/document/d/1aggwBCGu1CfsduOoKI0e6TRYOGtwwSceSBdKKkjCisA/edit?usp=sharing)
tutorials are available, covering:

- Kaptive's features and their scientific rationale
- How to run Kaptive
- Examples, illustrating how to run and interpret results
- Further investigations (e.g. exploring novel loci, IS insertions)

!!! note
    The tutorials are based on Kaptive 2.0, but the principles are similar
    for Kaptive 3.0.


## People

- [Dr. Tom Stanton](https://wyreslab.com/)
- [A/Prof. Kelly L. Wyres](https://wyreslab.com/)
- [Prof. Kathryn E. Holt](https://holtlab.net)
- [Dr. Ryan R. Wick](https://rrwick.github.io/)
- [A/Prof. Johanna
  Kenyon](https://experts.griffith.edu.au/45350-johanna-kenyon)

[Contact Kelly and Tom](mailto:kaptive.typing@gmail.com) for help with Kaptive,
or to report bugs or request features.

## References

[^1]: Van Der Graaf-van Bloois, L., Chen, H., Wagenaar, J.A., Zomer,
    A.L., 2023. Development of Kaptive databases for *Vibrio
    parahaemolyticus* O- and K-antigen genotyping. *Microbial Genomics*
    9. <https://doi.org/10.1099/mgen.0.001007>

[^2]: Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and
    Accurate in silico Antigen Typing with Kaptive 3.
    2025 _Microbial Genomics_ 11(6):001428.
    <https://doi.org/10.1099/mgen.0.001428>

[^3]: Wyres, K.L., Wick, R.R., Gorrie, C., Jenney, A., Follador, R.,
    Thomson, N.R., Holt, K.E., 2016. Identification of *Klebsiella*
    capsule synthesis loci from whole genome data. *Microbial Genomics*
    2, e000102. <https://doi.org/10.1099/mgen.0.000102>

[^4]: Lam, M.M.C., Wick, R.R., Judd, L.M., Holt, K.E., Wyres, K.L.,
    2022. Kaptive 2.0: updated capsule and lipopolysaccharide locus
    typing for the *Klebsiella pneumoniae* species complex. *Microbial
    Genomics* 8. <https://doi.org/10.1099/mgen.0.000800>

[^5]: Wick, R.R., Heinz, E., Holt, K.E., Wyres, K.L., 2018. Kaptive Web:
    User-Friendly Capsule and Lipopolysaccharide Serotype Prediction for
    *Klebsiella* Genomes. *Journal of Clinical Microbiology* 56,
    10.1128/jcm.00197-18. <https://doi.org/10.1128/jcm.00197-18>

[^6]: Wyres, K.L., Cahill, S.M., Holt, K.E., Hall, R.M., Kenyon, J.J.,
    2020. Identification of *Acinetobacter baumannii* loci for capsular
    polysaccharide (KL) and lipooligosaccharide outer core (OCL)
    synthesis in genome assemblies using curated reference databases
    compatible with Kaptive. *Microbial Genomics* 6.
    <https://doi.org/10.1099/mgen.0.000339>

[^7]: Cahill, S.M., Hall, R.M., Kenyon, J.J., 2022. An update to the
    database for *Acinetobacter baumannii* capsular polysaccharide locus
    typing extends the extensive and diverse repertoire of genes found
    at and outside the K locus. *Microbial Genomics* 8.
    <https://doi.org/10.1099/mgen.0.000878>
