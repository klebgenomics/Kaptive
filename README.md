<div align="center">
  <img src="docs/assets/logo.png" alt="Kaptive" width="400">
  <p>The tool for <em>in silico</em> serotyping</p>
  <br></br>
  <a href="https://pypi.org/project/kaptive/">
    <img src="https://img.shields.io/pypi/pyversions/kaptive" alt="PyPI - Python Version">
  </a>&nbsp;
  <a href="https://anaconda.org/bioconda/kaptive">
    <img src="https://img.shields.io/conda/vn/bioconda/kaptive" alt="Conda Version">
  </a>&nbsp;
  <a href="https://github.com/astral-sh/ruff">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="ruff">
  </a>&nbsp;
  <a href="https://github.com/astral-sh/ty">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ty/main/assets/badge/v0.json" alt="ty">
  </a>&nbsp;
</div>

For full **documentation**, including install and usage instructions, **click [here](https://klebgenomics.github.io/Kaptive)**.

## :arrow_forward: Tutorial

Step-by-step [video](https://klebnet.org/training/) and
[documented](https://docs.google.com/document/d/1aggwBCGu1CfsduOoKI0e6TRYOGtwwSceSBdKKkjCisA/edit?usp=sharing)
tutorials are available, covering:

- Kaptive's features and their scientific rationale
- How to run Kaptive
- Examples, illustrating how to run and interpret results
- Further investigations (e.g. exploring novel loci, IS insertions)

> **Note**: The tutorials are based on Kaptive 2.0, but the principles are similar
> for Kaptive 3.0.

## :mortar_board: Citation

If you use **Kaptive** in your work, please cite:

```bibtex
@article{mbs:/content/journal/mgen/10.1099/mgen.0.001428,
   author = "Stanton, Thomas David and Hetland, Marit A.K. and Löhr, Iren H. and Holt, Kathryn E. and Wyres, Kelly L.",
   title = "Fast and accurate in silico antigen typing with Kaptive 3", 
   journal= "Microbial Genomics",
   year = "2025",
   volume = "11",
   number = "6",
   pages = "",
   doi = "https://doi.org/10.1099/mgen.0.001428",
   url = "https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001428",
   publisher = "Microbiology Society",
   issn = "2057-5858",
   type = "Journal Article",
   keywords = "serotyping",
   keywords = "Klebsiella",
   keywords = "tools",
   keywords = "capsule",
   keywords = "sero-epidemiology",
   keywords = "antigen",
    eid = "001428",
   abstract = "Surface polysaccharides are common antigens in priority pathogens and therefore attractive targets for novel control strategies such as vaccines, monoclonal antibody and phage therapies. Distinct serotypes correspond to diverse polysaccharide structures that are encoded by distinct biosynthesis gene clusters; e.g. the Klebsiella pneumoniae species complex (KpSC) K- and O-loci encode the synthesis machinery for the capsule (K) and outer-lipopolysaccharides (O), respectively. We previously presented Kaptive and Kaptive 2, programmes to identify K- and O-loci directly from KpSC genome assemblies (later adapted for Acinetobacter baumannii), enabling sero-epidemiological analyses to guide vaccine and phage therapy development. However, for some KpSC genome collections, Kaptive (v≤2) was unable to type a high proportion of K-loci. Here, we identify the cause of this issue as assembly fragmentation and present a new version of Kaptive (v3) to circumvent this problem, reduce processing times and simplify output interpretation. We compared the performance of Kaptive v2 and Kaptive v3 for typing genome assemblies generated from subsampled Illumina read sets (decrements of 10× depth), for which a corresponding high-quality completed genome was also available to determine the ‘true’ loci (n=549 KpSC, n=198 A. baumannii). Both versions of Kaptive showed high rates of agreement to the matched true locus amongst ‘typeable’ locus calls (≥96% for ≥20× read depth), but Kaptive v3 was more sensitive, particularly for low-depth assemblies (at &lt;40× depth, v3 ranged 0.85–1 vs v2 0.09–0.94) and/or typing KpSC K-loci (e.g. 0.97 vs 0.82 for non-subsampled assemblies). Overall, Kaptive v3 was also associated with a higher rate of optimal outcomes; i.e. loci matching those in the reference database were correctly typed, and genuine novel loci were reported as untypeable (73–98% for v3 vs 7–77% for v2 for KpSC K-loci). Kaptive v3 was &gt;1 order of magnitude faster than Kaptive v2, making it easy to analyse thousands of assemblies on a desktop computer, facilitating broadly accessible in silico serotyping that is both accurate and sensitive. The Kaptive v3 source code is freely available on GitHub (https://github.com/klebgenomics/Kaptive), and has been implemented in Kaptive Web (https://kaptive-web.erc.monash.edu/).",
  }
```
