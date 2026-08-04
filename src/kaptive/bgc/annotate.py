r"""Gene prediction and functional annotation for Biosynthetic Gene Clusters (BGCs).

This module provides components for predicting Open Reading Frames (ORFs) from genome
assemblies using PyFGS, indexing translated protein sequences via strobemers, searching
against reference protein databases, and formatting annotation results into BED files.

Classes:
    Genes: Structure-of-Arrays (SoA) container storing interval coordinates, translated
        protein sequences, and contig assignments for predicted genes
        ([`Genes`][kaptive.bgc.annotate.Genes]).
    AnnotationResult: Container storing annotation outputs including gene SoA container,
        Randstrobe index, match seeds, database hit masks, top hit names/scores, and BED export
        ([`AnnotationResult`][kaptive.bgc.annotate.AnnotationResult]).
    Annotator: Parallel execution engine for ORF prediction, strobemer indexing, and database
        similarity searching ([`Annotator`][kaptive.bgc.annotate.Annotator]).
"""

from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

import numpy as np
import numpy.typing as npt
import pyfgs

from kaptive.core.collections import BatchedContainer
from kaptive.core.genome import GenomeAssembly
from kaptive.core.interval import Intervals
from kaptive.core.kmers import RandstrobeIndex, Seeds
from kaptive.core.pairwise import PairwiseAligner, PairwiseAlignments
from kaptive.core.seq import Sequences
from kaptive.db import Database


@dataclass(slots=True, frozen=True)
class Genes(BatchedContainer[Any, "Genes"]):
    r"""Structure-of-Arrays (SoA) container for predicted genes.

    Stores interval coordinates, translated amino acid sequences, and contig index
    mappings for a collection of predicted Open Reading Frames (ORFs). Inherits
    vectorised slicing and concatenation interface from
    [`BatchedContainer`][kaptive.core.collections.BatchedContainer].

    Attributes:
        intervals (Intervals): Genomic interval coordinates for predicted ORFs
            ([`Intervals`][kaptive.core.interval.Intervals]).
        translations (Sequences): Translated amino acid sequences
            ([`Sequences`][kaptive.core.seq.Sequences]).
        contig_indices (npt.NDArray[np.uint32]): 1D array of 0-based contig index
            offsets corresponding to parent contig names.
    """

    intervals: Intervals
    translations: Sequences
    contig_indices: npt.NDArray[np.uint32]

    def __len__(self) -> int:
        r"""Return total number of predicted genes in container.

        Returns:
            int: Number of genes.
        """
        return len(self.intervals)

    def __getitem__(self, item: int | slice | npt.NDArray[Any] | list[int]) -> Any:  # noqa: ANN401
        r"""Retrieve gene entry by index or slice container by range or boolean mask.

        Args:
            item (int | slice | npt.NDArray | list): Integer index, slice, boolean mask
                array, or index list.

        Returns:
            Any: Tuple of `(Interval, SeqRecord, int)` for integer index, or a new
                [`Genes`][kaptive.bgc.annotate.Genes] container for slices or arrays.
        """
        if isinstance(item, (int, np.integer)):
            return (
                self.intervals[item],
                self.translations[item],
                self.contig_indices[item],
            )
        return Genes(
            intervals=self.intervals[item],  # type: ignore
            translations=self.translations[item],  # type: ignore
            contig_indices=self.contig_indices[item],
        )

    @classmethod
    def empty(cls) -> "Genes":
        r"""Create an empty [`Genes`][kaptive.bgc.annotate.Genes] container.

        Returns:
            Genes: An empty gene container with zero elements.
        """
        return cls(
            intervals=Intervals.empty(),
            translations=Sequences.empty(),
            contig_indices=np.empty(0, dtype=np.uint32),
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenate multiple [`Genes`][kaptive.bgc.annotate.Genes] batches into a single container.

        Args:
            batches (Iterable[Genes]): Iterable of gene containers to concatenate.

        Returns:
            Genes: Combined gene container containing concatenated entries from all batches.
        """
        batches_list = list(batches)
        if not batches_list:
            return cls.empty()  # type: ignore
        return cls(
            intervals=Intervals.concat([b.intervals for b in batches_list]),
            translations=Sequences.concat([b.translations for b in batches_list]),
            contig_indices=np.concatenate([b.contig_indices for b in batches_list]),
        )


@dataclass(slots=True)
class AnnotationResult:
    r"""Container storing the results of genome annotation against a reference database.

    Holds predicted genes ([`Genes`][kaptive.bgc.annotate.Genes]), query protein strobemer
    index ([`RandstrobeIndex`][kaptive.core.kmers.RandstrobeIndex]), database hit seeds
    ([`Seeds`][kaptive.core.kmers.Seeds]), boolean hit mask, top hit metadata, and optional
    pairwise alignment results ([`PairwiseAlignments`][kaptive.core.pairwise.PairwiseAlignments]).

    Attributes:
        genes (Genes): SoA container of predicted genes
            ([`Genes`][kaptive.bgc.annotate.Genes]).
        translations_idx (RandstrobeIndex): Strobemer index constructed over translated gene
            sequences ([`RandstrobeIndex`][kaptive.core.kmers.RandstrobeIndex]).
        seeds (Seeds): Database hit seed matches ([`Seeds`][kaptive.core.kmers.Seeds]).
        hits_mask (npt.NDArray[np.bool_]): Boolean array indicating which genes matched database entries.
        top_hit_names (npt.NDArray[np.object_]): Array of top matching reference protein identifiers for each gene.
        top_hit_scores (npt.NDArray[np.float32]): Array of alignment/similarity scores for top database hits.
        contig_names (tuple[str, ...]): Tuple of contig names matching contig index offsets in genes.
        alignments (PairwiseAlignments | None): Optional pairwise alignment results
            ([`PairwiseAlignments`][kaptive.core.pairwise.PairwiseAlignments]), or `None` if disabled.
    """

    genes: Genes
    translations_idx: RandstrobeIndex
    seeds: Seeds
    hits_mask: npt.NDArray[np.bool_]
    top_hit_names: npt.NDArray[np.object_]  # Array of strings for each gene
    top_hit_scores: npt.NDArray[np.float32]  # Array of scores for each gene
    contig_names: tuple[str, ...]
    alignments: PairwiseAlignments | None = None

    def write_bed(self, path: str | Path, hits_only: bool = True) -> None:
        r"""Write predicted gene annotations to an extended 7-column BED file.

        Exports gene interval coordinates, strand orientation, and database hit metadata
        (top hit name and score in a 7th key=value column format) to a target file path.

        Args:
            path (str | Path): Destination file path for BED output.
            hits_only (bool): If `True`, only output genes with positive database hits. Default `True`.

        Raises:
            OSError: If writing to `path` fails due to filesystem or permissions errors.
        """
        mask = self.hits_mask if hits_only else np.ones(len(self.genes), dtype=bool)

        valid_indices = np.where(mask)[0]
        if len(valid_indices) == 0:
            with open(path, "w") as f:
                pass
            return

        filtered_genes = self.genes[valid_indices]
        c_names = [self.contig_names[i] for i in filtered_genes.contig_indices]

        starts = filtered_genes.intervals.starts
        ends = filtered_genes.intervals.ends
        strands = filtered_genes.intervals.strands

        strand_map = {1: "+", -1: "-", 0: "."}

        with open(path, "w") as f:
            for i, global_idx in enumerate(valid_indices):
                c_name = c_names[i]
                start = starts[i]
                end = ends[i]
                strand_char = strand_map.get(int(strands[i]), ".")
                name = str(global_idx)

                # Format tags for 7th column if there's a hit
                if self.hits_mask[global_idx]:
                    top_name = self.top_hit_names[global_idx]
                    top_score = self.top_hit_scores[global_idx]
                    tags = f"top_hit={top_name};score={top_score:.2f}"
                else:
                    tags = "."

                f.write(f"{c_name}\t{start}\t{end}\t{name}\t0\t{strand_char}\t{tags}\n")


class Annotator:
    r"""Parallel gene prediction and reference database annotation engine.

    Coordinates parallel ORF prediction on contigs using PyFGS, constructs query strobemer
    indexes ([`RandstrobeIndex`][kaptive.core.kmers.RandstrobeIndex]), queries a reference
    [`Database`][kaptive.db.Database] for matching hits, and optionally refines alignment
    scores using [`PairwiseAligner`][kaptive.core.pairwise.PairwiseAligner].

    Attributes:
        align (bool): Flag enabling optional pairwise sequence alignment for hits.
        aligner (PairwiseAligner | None): Pairwise alignment engine
            ([`PairwiseAligner`][kaptive.core.pairwise.PairwiseAligner]) if `align` is `True`, else `None`.
        whole_genome (bool): Flag configuring PyFGS gene finder mode (complete whole genome vs fragment mode).
    """

    def __init__(
        self,
        db: Database,
        align: bool = False,
        aligner_kwargs: dict[str, Any] | None = None,
        whole_genome: bool = False,
    ) -> None:
        r"""Initialize an [`Annotator`][kaptive.bgc.annotate.Annotator] instance.

        Args:
            db (Database): Target reference database ([`Database`][kaptive.db.Database])
                containing reference protein sequences.
            align (bool): If `True`, perform dynamic programming alignment on match seeds.
                Default `False`.
            aligner_kwargs (dict | None): Keyword arguments passed to
                [`PairwiseAligner`][kaptive.core.pairwise.PairwiseAligner]. Default `None`.
            whole_genome (bool): If `True`, treat input contigs as complete whole genomes in
                PyFGS gene finder. Default `False`.
        """
        self._db = db
        self.align = align
        self.aligner = PairwiseAligner(**(aligner_kwargs or {})) if align else None
        self.whole_genome = whole_genome

        # Build index for the database proteins using Randstrobe (must be sorted by hash for top_hits binary search)
        self._db_idx = RandstrobeIndex.build(self._db.translations, sort_by_hash=True)

    def __call__(self, genome: GenomeAssembly) -> AnnotationResult:
        r"""Run gene prediction and database annotation on a genome assembly.

        Predicts ORFs across all contigs in parallel using PyFGS, extracts protein
        translations, constructs strobemer indexes, queries database top hits, and returns
        an [`AnnotationResult`][kaptive.bgc.annotate.AnnotationResult].

        Args:
            genome (GenomeAssembly): Target input genome assembly
                ([`GenomeAssembly`][kaptive.core.genome.GenomeAssembly]).

        Returns:
            AnnotationResult: Complete annotation output container
                ([`AnnotationResult`][kaptive.bgc.annotate.AnnotationResult]).
        """
        model = pyfgs.Model.Complete
        finder = pyfgs.GeneFinder(model, whole_genome=self.whole_genome)

        offsets = genome.contigs.offsets
        lengths = genome.contigs.lengths
        seqs = genome.contigs.seqs

        # Collect sequence bytes for contigs
        seq_bytes_list = [seqs[offsets[i] : offsets[i] + lengths[i]].tobytes() for i in range(len(genome.contigs))]

        def _predict(seq_bytes: bytes) -> list[pyfgs.Gene]:
            if len(seq_bytes) < 3:
                return []
            return finder.find_genes(seq_bytes)

        # Predict ORFs across all contigs in parallel using ThreadPoolExecutor
        with ThreadPoolExecutor() as executor:
            batch_results = list(executor.map(_predict, seq_bytes_list))

        starts: list[int] = []
        ends: list[int] = []
        strands: list[int] = []
        contig_indices: list[int] = []
        protein_seqs: list[bytes] = []

        for contig_idx, genes in enumerate(batch_results):
            for g in genes:
                starts.append(g.start)
                ends.append(g.end)
                strands.append(g.strand)
                contig_indices.append(contig_idx)
                protein_seqs.append(g.translation())

        intervals = Intervals(
            starts=np.array(starts, dtype=np.int32),
            ends=np.array(ends, dtype=np.int32),
            strands=np.array(strands, dtype=np.int8),
        )

        genes_soa = Genes(
            intervals=intervals,
            translations=Sequences.from_bytes(protein_seqs),
            contig_indices=np.array(contig_indices, dtype=np.uint32),
        )

        # Build query translation index
        translations_idx = RandstrobeIndex.build(genes_soa.translations)

        # Compare translations against database to find seeds
        seeds = self._db_idx.top_hits(translations_idx)

        # Build hits mask and track top hit info
        num_genes = len(genes_soa)
        hits_mask = np.zeros(num_genes, dtype=bool)
        top_hit_names = np.full(num_genes, "", dtype=object)
        top_hit_scores = np.zeros(num_genes, dtype=np.float32)

        if len(seeds) > 0:
            hits_mask[seeds.query_indices] = True

            for q_idx, t_idx, score in zip(seeds.query_indices, seeds.target_indices, seeds.scores):
                top_hit_names[q_idx] = self._db.translations.ids[t_idx]
                top_hit_scores[q_idx] = float(score)

        alignments = None
        if self.aligner is not None and len(seeds) > 0:
            q_seqs, t_seqs = seeds.extract_sequences(genes_soa.translations, self._db.translations)
            alignments = self.aligner(q_seqs, t_seqs, seeds=seeds)

            # Update scores if aligned
            for i, q_idx in enumerate(seeds.query_indices):
                top_hit_scores[q_idx] = float(alignments.scores[i])

        return AnnotationResult(
            genes=genes_soa,
            translations_idx=translations_idx,
            seeds=seeds,
            hits_mask=hits_mask,
            top_hit_names=top_hit_names,
            top_hit_scores=top_hit_scores,
            contig_names=genome.contigs.ids,
            alignments=alignments,
        )
