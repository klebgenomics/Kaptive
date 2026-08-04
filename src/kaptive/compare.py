r"""Multi-locus comparative genomic analysis and sequence alignment containers.

This module provides high-performance Structure-of-Arrays (SoA) data containers
and an indexed comparison engine for performing all-vs-all forward protein alignments
across multiple genomic loci.

Classes:
    [`LocusComparisonEdges`][kaptive.compare.LocusComparisonEdges]: Container storing
        pairwise cross-locus protein alignments.
    [`LocusComparisons`][kaptive.compare.LocusComparisons]: Complete result set for multi-locus comparison.
    [`LocusData`][kaptive.compare.LocusData]: Generalised input data container representing a single locus.
    [`LocusComparator`][kaptive.compare.LocusComparator]: Vectorised engine for multi-locus pairwise comparisons.
"""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Self

import numpy as np
import numpy.typing as npt

from kaptive.core.collections import BatchedContainer
from kaptive.core.interval import Intervals
from kaptive.core.kmers import RandstrobeIndex
from kaptive.core.pairwise import PairwiseAligner, PairwiseAlignments
from kaptive.core.seq import Sequences

if TYPE_CHECKING:
    from kaptive.serotyping import LocusPieces


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class LocusComparisonEdges(BatchedContainer[Any, "LocusComparisonEdges"]):
    r"""A high-performance SoA container for forward cross-locus protein alignments.

    Stores index pointers and alignment statistics for protein-level pairwise hits
    identified between different genomic loci. Inherits batch slicing and concatenation
    capabilities from [`BatchedContainer`][kaptive.core.collections.BatchedContainer].

    Attributes:
        query_locus_indices (npt.NDArray[np.int32]): Locus-level indices for query loci.
        target_locus_indices (npt.NDArray[np.int32]): Locus-level indices for target loci.
        query_indices (npt.NDArray[np.int32]): Protein indices within the query locus.
        target_indices (npt.NDArray[np.int32]): Protein indices within the target locus.
        global_query_indices (npt.NDArray[np.int32]): Global protein indices across concatenated loci.
        global_target_indices (npt.NDArray[np.int32]): Global protein indices across concatenated loci.
        alignments (PairwiseAlignments): Pairwise alignment metrics for each edge.
    """

    # 1. Pointers to the locus-level sequences (useful for metadata)
    query_locus_indices: npt.NDArray[np.int32]
    target_locus_indices: npt.NDArray[np.int32]

    # 2. Pointers to the protein-level sequences within each locus
    query_indices: npt.NDArray[np.int32]
    target_indices: npt.NDArray[np.int32]

    # 3. Global pointers (if all loci were concatenated)
    global_query_indices: npt.NDArray[np.int32]
    global_target_indices: npt.NDArray[np.int32]

    # 4. The actual alignment results for these edges
    alignments: PairwiseAlignments

    def __len__(self) -> int:
        r"""Return the total number of cross-locus alignment edges.

        Returns:
            int: The number of alignment edges stored in the container.
        """
        return len(self.query_locus_indices)

    def __getitem__(self, item: int | slice | npt.NDArray[Any] | list[int]) -> "Any | LocusComparisonEdges":
        r"""Slice or filter alignment edges using a slice or array mask.

        Args:
            item (int | slice | npt.NDArray | list): Slice or boolean/integer mask array.

        Returns:
            Any | LocusComparisonEdges: A new container with sliced alignment edges.

        Raises:
            NotImplementedError: If `item` is a single integer.
        """
        if isinstance(item, (int, np.integer)):
            raise NotImplementedError("Single item access not implemented for LocusComparisonEdges")
        return LocusComparisonEdges(
            query_locus_indices=self.query_locus_indices[item],
            target_locus_indices=self.target_locus_indices[item],
            query_indices=self.query_indices[item],
            target_indices=self.target_indices[item],
            global_query_indices=self.global_query_indices[item],
            global_target_indices=self.global_target_indices[item],
            alignments=self.alignments[item],  # type: ignore
        )

    @classmethod
    def empty(cls) -> "LocusComparisonEdges":
        r"""Create an empty `LocusComparisonEdges` instance.

        Returns:
            LocusComparisonEdges: An empty edges container with zero-length arrays.
        """
        return cls(
            query_locus_indices=np.empty(0, dtype=np.int32),
            target_locus_indices=np.empty(0, dtype=np.int32),
            query_indices=np.empty(0, dtype=np.int32),
            target_indices=np.empty(0, dtype=np.int32),
            global_query_indices=np.empty(0, dtype=np.int32),
            global_target_indices=np.empty(0, dtype=np.int32),
            alignments=PairwiseAlignments.empty(),
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenate multiple `LocusComparisonEdges` batches into a single container.

        Args:
            batches (list[LocusComparisonEdges]): List of edge containers to merge.

        Returns:
            LocusComparisonEdges: Concatenated container holding all input batches.
        """
        if not batches:
            return cls.empty()  # type: ignore
        return cls(
            query_locus_indices=np.concatenate([b.query_locus_indices for b in batches]),
            target_locus_indices=np.concatenate([b.target_locus_indices for b in batches]),
            query_indices=np.concatenate([b.query_indices for b in batches]),
            target_indices=np.concatenate([b.target_indices for b in batches]),
            global_query_indices=np.concatenate([b.global_query_indices for b in batches]),
            global_target_indices=np.concatenate([b.global_target_indices for b in batches]),
            alignments=PairwiseAlignments.concat([b.alignments for b in batches]),
        )


@dataclass(slots=True, frozen=True)
class LocusComparisons:
    r"""Complete result container for multi-locus comparative analysis.

    Holds pairwise alignment edges between loci, global gene metadata, locus offsets,
    and physical genomic intervals normalized for visualization.

    Attributes:
        edges (LocusComparisonEdges): Alignment edges between compared loci.
        locus_names (tuple[str, ...]): Names of the compared loci.
        locus_lengths (npt.NDArray[np.int32]): Protein counts for each locus.
        locus_offsets (npt.NDArray[np.int32]): Global gene index offsets for each locus.
        gene_names (npt.NDArray[np.object_]): Flattened array of gene identifiers.
        gene_descriptions (npt.NDArray[np.object_]): Flattened array of gene descriptions.
        gene_states (npt.NDArray[np.int8]): Flattened array of gene classification states.
        gene_intervals (Intervals): Pre-normalised physical gene intervals for plotting.
    """

    edges: LocusComparisonEdges

    # Metadata for the loci being compared
    locus_names: tuple[str, ...]
    locus_lengths: npt.NDArray[np.int32]
    locus_offsets: npt.NDArray[np.int32]

    # Global flattened arrays for all genes across all loci
    gene_names: npt.NDArray[np.object_]
    gene_descriptions: npt.NDArray[np.object_]
    gene_states: npt.NDArray[np.int8]

    # The physical intervals of all genes, pre-normalised for plotting
    gene_intervals: Intervals


@dataclass(slots=True, frozen=True)
class LocusData:
    r"""A generalised container representing a single locus for comparison.

    Attributes:
        proteins (Sequences): Protein sequences within the locus.
        name (str): Locus identifier or name.
        backbone (Intervals): Physical genomic coordinates of genes in the locus.
        pieces (LocusPieces | None): Segmented contig piece layout, if fragmented.
        gene_ctg_indices (npt.NDArray[np.uint32] | None): Contig assignment indices for genes.
        gene_states (npt.NDArray[np.int8] | None): Array of GeneState integer values for genes in locus.
        gene_descriptions (npt.NDArray[np.object_] | Sequence[str] | None): Product description strings for genes.
    """

    proteins: Sequences
    name: str
    backbone: Intervals
    pieces: "LocusPieces | None" = None
    gene_ctg_indices: npt.NDArray[np.uint32] | None = None
    gene_states: npt.NDArray[np.int8] | None = None
    gene_descriptions: npt.NDArray[np.object_] | Sequence[str] | None = None


class LocusComparator:
    r"""Vectorised engine for multi-locus all-vs-all forward protein comparison.

    Uses strobemer indexing ([`RandstrobeIndex`][kaptive.core.kmers.RandstrobeIndex]) to quickly identify
    candidate homology hits between loci, and aligns candidates using
    [`PairwiseAligner`][kaptive.core.pairwise.PairwiseAligner].
    """

    def __init__(
        self,
        k: int = 10,
        s: int = 5,
        min_score: int = 1,
        aligner_kwargs: dict | None = None,  # type: ignore
    ) -> None:
        r"""Initialise the multi-locus comparator engine.

        Args:
            k (int): Strobemer seed length. Defaults to 10.
            s (int): Strobemer sampling stride. Defaults to 5.
            min_score (int): Minimum seed hit score threshold. Defaults to 1.
            aligner_kwargs (dict | None): Optional arguments passed to
                [`PairwiseAligner`][kaptive.core.pairwise.PairwiseAligner].
        """
        self.k = k
        self.s = s
        self.min_score = min_score
        self.aligner = PairwiseAligner(**(aligner_kwargs or {}))

    def __call__(self, inputs: Sequence[LocusData]) -> LocusComparisons:
        r"""Perform multi-locus pairwise comparisons across all input loci.

        Args:
            inputs (Sequence[LocusData]): Sequence of input loci data containers to compare.

        Returns:
            LocusComparisons: Complete comparison results including pairwise edges and normalized intervals.
        """
        loci = [inp.proteins for inp in inputs]
        locus_names = [inp.name for inp in inputs]
        backbones = [inp.backbone for inp in inputs]
        locus_pieces = [inp.pieces for inp in inputs]
        gene_ctg_indices = [inp.gene_ctg_indices for inp in inputs]

        n_loci = len(loci)

        # Build global sequences to extract metadata
        global_seqs = Sequences.concat(loci) if n_loci > 0 else Sequences.empty()
        gene_names = np.array(global_seqs.ids, dtype=object)

        desc_list = []
        state_list = []
        for inp in inputs:
            n_genes = len(inp.proteins)
            if len(inp.backbone) != n_genes:
                raise ValueError(
                    f"Locus '{inp.name}': backbone length ({len(inp.backbone)}) does not match protein count ({n_genes})"
                )

            if inp.gene_descriptions is not None:
                raw_desc = np.asarray(inp.gene_descriptions)
                if raw_desc.dtype.kind in ("S", "a"):
                    decoded_desc = np.char.decode(raw_desc, "utf-8")
                    d_arr = np.asarray(decoded_desc, dtype=object)
                elif raw_desc.dtype == object or any(isinstance(x, (bytes, np.bytes_)) for x in raw_desc.flat):
                    decoded_list = [
                        x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x) if x is not None else ""
                        for x in raw_desc.flat
                    ]
                    d_arr = np.asarray(decoded_list, dtype=object).reshape(raw_desc.shape)
                else:
                    d_arr = np.asarray(raw_desc, dtype=object)

                if len(d_arr) != n_genes:
                    raise ValueError(
                        f"Locus '{inp.name}': gene_descriptions length ({len(d_arr)}) does not match protein count ({n_genes})"
                    )
                desc_list.append(d_arr)
            else:
                desc_list.append(np.array([""] * n_genes, dtype=object))

            if inp.gene_states is not None:
                s_arr = np.asarray(inp.gene_states, dtype=np.int8)
                if len(s_arr) != n_genes:
                    raise ValueError(
                        f"Locus '{inp.name}': gene_states length ({len(s_arr)}) does not match protein count ({n_genes})"
                    )
                state_list.append(s_arr)
            else:
                from kaptive.serotyping.models import GeneState

                state_list.append(np.full(n_genes, GeneState.NORMAL.value, dtype=np.int8))

        if n_loci > 0:
            gene_descriptions = np.concatenate(desc_list, dtype=object)
            gene_states = np.concatenate(state_list, dtype=np.int8)
        else:
            gene_descriptions = np.empty(0, dtype=object)
            gene_states = np.empty(0, dtype=np.int8)

        # Normalise backbones
        norm_backbones = []
        for i, bb in enumerate(backbones):
            if locus_pieces is not None and i < len(locus_pieces) and locus_pieces[i] is not None:
                lp = locus_pieces[i]
                # Determine piece indices for each gene
                p_idx = np.zeros(len(bb), dtype=np.int32)
                for p in range(len(lp)):  # type: ignore
                    mask = (bb.starts >= lp.starts[p]) & (bb.ends <= lp.ends[p])  # type: ignore
                    if gene_ctg_indices is not None and gene_ctg_indices[i] is not None:
                        mask &= gene_ctg_indices[i] == lp.ctg_indices[p]  # type: ignore
                    p_idx[mask] = p
                # Assume piece_order is linear 0..N-1 as Serotyper.__call__ builds them in expected order
                p_order = np.arange(len(lp), dtype=np.int32)  # type: ignore
                bb_norm = bb.arrange(p_idx, p_order, lp.starts, lp.ends, lp.strands)  # type: ignore
                norm_backbones.append(bb_norm)
            else:
                if len(bb) > 0:
                    bb_norm = bb.shift(-np.min(bb.starts))  # type: ignore
                else:
                    bb_norm = bb
                norm_backbones.append(bb_norm)

        # Construct global intervals
        if norm_backbones:
            global_intervals = Intervals(
                starts=np.concatenate([b.starts for b in norm_backbones]),
                ends=np.concatenate([b.ends for b in norm_backbones]),
                strands=np.concatenate([b.strands for b in norm_backbones]),
                original_indices=np.concatenate([b.original_indices for b in norm_backbones]),
            )
        else:
            global_intervals = Intervals(
                np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int8)
            )

        # We also want to track global offsets for the proteins
        locus_lengths = np.array([len(seq_locus) for seq_locus in loci], dtype=np.int32)
        locus_offsets = np.zeros(n_loci, dtype=np.int32)
        if n_loci > 1:
            np.cumsum(locus_lengths[:-1], out=locus_offsets[1:])

        if n_loci <= 1:
            edges = LocusComparisonEdges.empty()
            return LocusComparisons(
                edges=edges,
                locus_names=tuple(locus_names),
                locus_lengths=locus_lengths,
                locus_offsets=locus_offsets,
                gene_names=gene_names,
                gene_descriptions=gene_descriptions,
                gene_states=gene_states,
                gene_intervals=global_intervals,
            )

        # Build Randstrobe indices for all loci
        # Targets must be sorted for binary search lookups
        target_indices = [RandstrobeIndex.build(seq_locus, k=self.k, s=self.s, sort_by_hash=True) for seq_locus in loci]
        # Queries must NOT be sorted so they can be streamed positionally
        query_indices = [RandstrobeIndex.build(seq_locus, k=self.k, s=self.s, sort_by_hash=False) for seq_locus in loci]

        edge_batches = []

        # Forward upper-triangle comparisons only (no self comparisons, no reverse)
        for i in range(n_loci):
            for j in range(i + 1, n_loci):
                # We align proteins from locus i to their best hit in locus j
                seeds = target_indices[j].top_hits(query_indices[i], min_score=self.min_score)
                if len(seeds) == 0:
                    continue

                alignments = self.aligner.align_seeds(loci[i], loci[j], seeds)

                n_edges = len(seeds)

                # Construct SoA pointers
                batch = LocusComparisonEdges(
                    query_locus_indices=np.full(n_edges, i, dtype=np.int32),
                    target_locus_indices=np.full(n_edges, j, dtype=np.int32),
                    query_indices=seeds.query_indices.astype(np.int32),
                    target_indices=seeds.target_indices.astype(np.int32),
                    global_query_indices=seeds.query_indices.astype(np.int32) + locus_offsets[i],
                    global_target_indices=seeds.target_indices.astype(np.int32) + locus_offsets[j],
                    alignments=alignments,
                )
                edge_batches.append(batch)

        if not edge_batches:
            edges = LocusComparisonEdges.empty()
        else:
            edges = LocusComparisonEdges.concat(edge_batches)

        return LocusComparisons(
            edges=edges,
            locus_names=tuple(locus_names),
            locus_lengths=locus_lengths,
            locus_offsets=locus_offsets,
            gene_names=gene_names,
            gene_descriptions=gene_descriptions,
            gene_states=gene_states,
            gene_intervals=global_intervals,
        )
