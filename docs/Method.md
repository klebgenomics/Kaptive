# Method

## kaptive assembly

For each input assembly, Kaptive runs the
`kaptive.assembly.typing_pipeline` which does the following:

1.  Aligns locus gene nucleotide sequences to the assembly contig
    sequences using [minimap2](https://github.com/lh3/minimap2).
2.  Identifies the best matching locus type using the [scoring
    algorithm](Method.md#scoring-algorithm).
3.  Extracts the locus gene sequences from the assembly contig
    sequences.
4.  Predicts the serotype/phenotype based on the gene content.

<a id="scoring-algorithm"></a>

### Scoring algorithm

1.  A matrix is initialised with one row per locus and one column per
    score metric.
2.  For each gene found, the best alignment is chosen and the scores are
    added to the respective locus row if the coverage passes the
    threshold (`--min-cov`).
3.  The matrix is then weighted by the `--weight-metric` and the scores
    are selected from the column corresponding to the `--score-metric`.
4.  The top N loci (`--n-best`) are selected to be fully aligned to the
    assembly.
5.  Steps 1 and 2 are repeated and the best locus is selected from the
    column corresponding to the `--score-metric`.

The score metric can be explicitly set by the flag `--score-metric`, the
options are:

- `0` (`AS`) - Sum of the alignment scores per query
- `1` (`mlen`) - Sum of the number of matching bases in the alignment
  per query
- `2` (`blen`) - Sum of the number of aligned bases per query
- `3` (`q_len`) - Sum of the number of bases of each query found
  (regardless of whether they are aligned)

The weighting can be explicitly set by the flag `--weight-metric`; the
options are:

- `0` - No weighting
- `1` - Number of genes found
- `2` - Number of genes expected
- `3` - Proportion of genes found
- `4` - Sum of the number of aligned bases per query (`blen`)
- `5` - Sum of the number of bases of each query found (regardless of
  whether they are aligned) (`q_len`)

!!! note
    The gene score matrix can be written to a TSV file using the
    `--scores` flag, however this will not type the assembly or
    reconstruct the locus.


<a id="locus-reconstruction"></a>

### Locus reconstruction

After the best matching locus type has been identified, Kaptive will:

1.  For each contig, the ranges from the full-length locus alignments of
    the best match are extracted.

2.  The ranges are merged together if they are within the distance of
    the largest locus in the database.

3.  The merged ranges are used to create `LocusPiece` objects and the
    sequence is extracted from the assembly contig.

4.  Gene alignments are culled twice to determine the gene content:  
    1.  The first removes alignments overlapping genes from the best
        match.
    2.  The remaining alignments that are not part of the best match are
        culled so that the best alignment is kept and each alignment
        represents the best gene for that part of the contig.

5.  For each remaining gene alignment, the gene is then evaluated.

!!! note
    You may see multiple entries for the same gene in the gene details
    columns. This can happen when there are multiple significant
    alignments for a gene that has been split, either over contigs or by
    an artefact such as and insertion element. This will be evident by the
    low percent coverage of the entries, which should all add up to ~100%
    if they are indeed part of the same gene.


<a id="gene-evaluation"></a>

### Gene evaluation

For each `GeneResult` object, Kaptive will:

1.  Check whether the gene is **partial** by determining if the gene
    overlaps the contig boundaries.
2.  Extract the DNA sequence from the assembly contig and translate to
    amino acid.
3.  Perform pairwise alignment to the reference gene amino acid sequence
    and calculate percent identity.
4.  Check for **truncation** by determining if the amino acid sequence
    length is \<95% of the reference gene protein length.
5.  Determine whether the gene belongs to the biosynthetic gene cluster
    (**inside locus**) or not (**outside locus**).

!!! note
    Partial genes are *not* considered for truncation. This prevents false
    positive truncation calls in fragmented assemblies which may otherwise
    have an impact on phenotype prediction.


<a id="phenotype-prediction"></a>

### Phenotype prediction

As of Kaptive 3, we have added the ability to predict the resulting
phenotype of the assembly. This is similar to how the *Type* was
reported in previous versions, but now includes the ability to predict
specific phenotypes based on known mutations/modifications in a given
set of locus genes as defined in the database [logic
file](Databases.md#phenotype-logic).

### Confidence score

Kaptive with finally calculate how confident it is in the prediction,
which is explained [here](Interpreting-the-results.md#confidence-score).
