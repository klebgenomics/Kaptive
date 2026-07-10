from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from kaptive.core.interval import Intervals
from kaptive.core.pairwise import PairwiseAligner, PairwiseAlignments, RandstrobeIndex
from kaptive.core.seq import Sequences


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class LocusComparisonEdges:
    """A high-performance SoA container for all forward cross-locus protein alignments."""
    
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
        return len(self.query_locus_indices)

    @classmethod
    def empty(cls) -> 'LocusComparisonEdges':
        return cls(
            query_locus_indices=np.empty(0, dtype=np.int32),
            target_locus_indices=np.empty(0, dtype=np.int32),
            query_indices=np.empty(0, dtype=np.int32),
            target_indices=np.empty(0, dtype=np.int32),
            global_query_indices=np.empty(0, dtype=np.int32),
            global_target_indices=np.empty(0, dtype=np.int32),
            alignments=PairwiseAlignments.empty()
        )

    @classmethod
    def concat(cls, batches: list['LocusComparisonEdges']) -> 'LocusComparisonEdges':
        if not batches:
            return cls.empty()
        return cls(
            query_locus_indices=np.concatenate([b.query_locus_indices for b in batches]),
            target_locus_indices=np.concatenate([b.target_locus_indices for b in batches]),
            query_indices=np.concatenate([b.query_indices for b in batches]),
            target_indices=np.concatenate([b.target_indices for b in batches]),
            global_query_indices=np.concatenate([b.global_query_indices for b in batches]),
            global_target_indices=np.concatenate([b.global_target_indices for b in batches]),
            alignments=PairwiseAlignments.concat([b.alignments for b in batches])
        )


@dataclass(slots=True, frozen=True)
class LocusComparisons:
    edges: LocusComparisonEdges
    
    # Metadata for the loci being compared
    locus_names: tuple[str, ...]
    locus_lengths: npt.NDArray[np.int32]
    locus_offsets: npt.NDArray[np.int32]
    
    # Global flattened arrays for all genes across all loci
    gene_names: npt.NDArray[np.object_]
    gene_descriptions: npt.NDArray[np.object_]
    
    # The physical intervals of all genes, pre-normalised for plotting
    gene_intervals: Intervals


class LocusComparator:
    """
    Highly vectorized engine for comparing multiple loci (collections of protein sequences)
    in an all-vs-all forward manner.
    """
    def __init__(self, k: int = 10, s: int = 5, min_score: int = 1, aligner_kwargs: dict | None = None):
        self.k = k
        self.s = s
        self.min_score = min_score
        self.aligner = PairwiseAligner(**(aligner_kwargs or {}))

    def __call__(self, 
                 loci: Sequence[Sequences], 
                 locus_names: Sequence[str],
                 backbones: Sequence[Intervals],
                 locus_pieces: Sequence['LocusPieces'] | None = None,
                 gene_ctg_indices: Sequence[npt.NDArray[np.uint32]] | None = None) -> LocusComparisons:
        n_loci = len(loci)
        
        # Build global sequences to extract metadata
        global_seqs = Sequences.concat(loci) if n_loci > 0 else Sequences.empty()
        gene_names = np.array(global_seqs.ids, dtype=object)
        gene_descriptions = np.array([""] * len(global_seqs.ids), dtype=object)
            
        # Normalise backbones
        norm_backbones = []
        for i, bb in enumerate(backbones):
            if locus_pieces is not None and i < len(locus_pieces) and locus_pieces[i] is not None:
                lp = locus_pieces[i]
                # Determine piece indices for each gene
                p_idx = np.zeros(len(bb), dtype=np.int32)
                for p in range(len(lp)):
                    mask = (bb.starts >= lp.starts[p]) & (bb.ends <= lp.ends[p])
                    if gene_ctg_indices is not None and gene_ctg_indices[i] is not None:
                        mask &= (gene_ctg_indices[i] == lp.ctg_indices[p])
                    p_idx[mask] = p
                # Assume piece_order is linear 0..N-1 as Serotyper.__call__ builds them in expected order
                p_order = np.arange(len(lp), dtype=np.int32)
                bb_norm = bb.arrange(p_idx, p_order, lp.starts, lp.ends, lp.strands)
                norm_backbones.append(bb_norm)
            else:
                if len(bb) > 0:
                    bb_norm = bb.shift(-np.min(bb.starts))
                else:
                    bb_norm = bb
                norm_backbones.append(bb_norm)
            
        # Construct global intervals
        if norm_backbones:
            global_intervals = Intervals(
                starts=np.concatenate([b.starts for b in norm_backbones]),
                ends=np.concatenate([b.ends for b in norm_backbones]),
                strands=np.concatenate([b.strands for b in norm_backbones]),
                original_indices=np.concatenate([b.original_indices for b in norm_backbones])
            )
        else:
            global_intervals = Intervals(np.empty(0, dtype=np.int32),
                                         np.empty(0, dtype=np.int32),
                                         np.empty(0, dtype=np.int8))
            
        # We also want to track global offsets for the proteins
        locus_lengths = np.array([len(l) for l in loci], dtype=np.int32)
        locus_offsets = np.zeros(n_loci, dtype=np.int32)
        if n_loci > 1:
            np.cumsum(locus_lengths[:-1], out=locus_offsets[1:])
            
        if n_loci <= 1:
            edges = LocusComparisonEdges.empty()
            return LocusComparisons(edges, tuple(locus_names), locus_lengths, locus_offsets, gene_names,
                                    gene_descriptions, global_intervals)

        # Build Randstrobe indices for all loci
        # Targets must be sorted for binary search lookups
        target_indices = [RandstrobeIndex.build(l, k=self.k, s=self.s, sort_by_hash=True) for l in loci]
        # Queries must NOT be sorted so they can be streamed positionally
        query_indices = [RandstrobeIndex.build(l, k=self.k, s=self.s, sort_by_hash=False) for l in loci]
        
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
                    alignments=alignments
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
            gene_intervals=global_intervals
        )
