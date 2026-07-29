r"""Core data structures, container protocols, and computational engines.

The `kaptive.core` sub-package provides high-performance data structures,
interval operations, k-mer indexing algorithms, and pairwise sequence alignment
routines underlying Kaptive's serotyping and BGC prediction functionality.

Modules:
    alignment: CIGAR parsing, alignment models, and batched alignment containers
        (e.g., [`Alignments`][kaptive.core.alignment.Alignments],
        [`Alignment`][kaptive.core.alignment.Alignment],
        [`Cigars`][kaptive.core.alignment.Cigars]).
    collections: Protocol base classes for Structure-of-Arrays containers
        (e.g., [`BatchedContainer`][kaptive.core.collections.BatchedContainer],
        [`RaggedArrayContainer`][kaptive.core.collections.RaggedArrayContainer]).
    genome: Genome assembly parsing and FASTA sequence iteration
        (e.g., [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly],
        [`FastaReader`][kaptive.core.genome.FastaReader]).
    interval: Genomic interval representations, clustering, and overlap resolution
        (e.g., [`Interval`][kaptive.core.interval.Interval],
        [`Intervals`][kaptive.core.interval.Intervals],
        [`Strand`][kaptive.core.interval.Strand]).
    kmers: MinHash and Randstrobe k-mer indexing and matching engines
        (e.g., [`RandstrobeIndex`][kaptive.core.kmers.RandstrobeIndex],
        [`FracMinHashIndex`][kaptive.core.kmers.FracMinHashIndex],
        [`Seeds`][kaptive.core.kmers.Seeds]).
    pairwise: Batched dynamic programming pairwise alignment engines
        (e.g., [`PairwiseAligner`][kaptive.core.pairwise.PairwiseAligner],
        [`PairwiseAlignments`][kaptive.core.pairwise.PairwiseAlignments]).
    seq: Sequence storage, vectorised translation, and sequence record metadata
        (e.g., [`Sequences`][kaptive.core.seq.Sequences],
        [`SeqRecord`][kaptive.core.seq.SeqRecord]).
"""
