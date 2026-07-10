"""
Module for performing pairwise protein alignments using a numba-powered banded Gotoh Smith-Waterman kernel.
Also contains a Randstrobe strobemer generation and jaccard implementation to act as a pre-filter.
"""
from dataclasses import dataclass
from functools import cache
from typing import Union

import numpy as np
import numpy.typing as npt
from numba import njit, prange

from kaptive.core.seq import Sequences

# Constants ------------------------------------------------------------------------------------------------------------
_INF = np.int32(-1_000_000_000)  # Safe negative infinity for int32 DP arithmetic


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class PairwiseAlignment:
    """A lightweight, immutable container for the results of a single pairwise sequence alignment.

    This class holds the summary statistics of an alignment, such as the score and the counts of
    matches, mismatches, and gaps. It is typically created when accessing a single element from a
    `PairwiseAlignments` collection.

    Attributes:
        score (int): The final alignment score, calculated using a substitution matrix (e.g., BLOSUM62)
            and gap penalties.
        matches (int): The total number of matching characters in the alignment.
        mismatches (int): The total number of mismatched characters.
        gaps (int): The total number of gaps (insertions or deletions).
    q_start (int): The 0-based start coordinate on the query sequence.
    q_end (int): The 0-based end coordinate on the query sequence.
    t_start (int): The 0-based start coordinate on the target sequence.
    t_end (int): The 0-based end coordinate on the target sequence.
    """
    score: int
    matches: int
    mismatches: int
    gaps: int
    q_start: int
    q_end: int
    t_start: int
    t_end: int

    @property
    def pident(self) -> float:
        """Calculates the percent identity of the aligned region.

        Percent identity is defined as `(matches / (matches + mismatches + gaps)) * 100`.

        Returns:
            float: The percent identity, or 0.0 if the total alignment length is zero.
        """
        total = self.matches + self.mismatches + self.gaps
        return (self.matches / total) * 100.0 if total > 0 else 0.0


@dataclass(frozen=True, slots=True)
class PairwiseAlignments:
    """A high-performance SoA container for the results of multiple pairwise alignments.

    This class stores alignment statistics in a Structure-of-Arrays (SoA) layout using NumPy arrays.
    This design enables highly efficient, vectorized calculations on large sets of alignment results,
    such as calculating percent identity for all alignments at once.

    Attributes:
        scores (npt.NDArray[np.int32]): A 1D array of alignment scores.
        matches (npt.NDArray[np.int32]): A 1D array of match counts.
        mismatches (npt.NDArray[np.int32]): A 1D array of mismatch counts.
        gaps (npt.NDArray[np.int32]): A 1D array of gap counts.
        q_starts (npt.NDArray[np.int32]): A 1D array of query start coordinates.
        q_ends (npt.NDArray[np.int32]): A 1D array of query end coordinates.
        t_starts (npt.NDArray[np.int32]): A 1D array of target start coordinates.
        t_ends (npt.NDArray[np.int32]): A 1D array of target end coordinates.
    """
    scores:     npt.NDArray[np.int32]
    matches:    npt.NDArray[np.int32]
    mismatches: npt.NDArray[np.int32]
    gaps:       npt.NDArray[np.int32]
    q_starts:   npt.NDArray[np.int32]
    q_ends:     npt.NDArray[np.int32]
    t_starts:   npt.NDArray[np.int32]
    t_ends:     npt.NDArray[np.int32]

    def __len__(self) -> int:
        """Returns the number of alignments in the batch."""
        return len(self.scores)

    def to_dict(self) -> dict:
        """Returns a dictionary containing the SoA arrays for orjson serialization."""
        return {
            'scores': self.scores,
            'matches': self.matches,
            'mismatches': self.mismatches,
            'gaps': self.gaps,
            'q_starts': self.q_starts,
            'q_ends': self.q_ends,
            't_starts': self.t_starts,
            't_ends': self.t_ends
        }

    @classmethod
    def from_dict(cls, d: dict) -> 'PairwiseAlignments':
        """Deserializes alignments from a dictionary."""
        return cls(
            np.array(d['scores'], dtype=np.int32),
            np.array(d['matches'], dtype=np.int32),
            np.array(d['mismatches'], dtype=np.int32),
            np.array(d['gaps'], dtype=np.int32),
            np.array(d['q_starts'], dtype=np.int32),
            np.array(d['q_ends'], dtype=np.int32),
            np.array(d['t_starts'], dtype=np.int32),
            np.array(d['t_ends'], dtype=np.int32)
        )

    def __getitem__(self, item) -> Union[PairwiseAlignment, 'PairwiseAlignments']:
        """Accesses alignment results by index, slice, or boolean mask.

        - If `item` is an integer, it returns a single `PairwiseAlignment` object.
        - If `item` is a slice or mask, it returns a new, smaller `PairwiseAlignments` collection.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return PairwiseAlignment(
                score=self.scores[item],
                matches=self.matches[item],
                mismatches=self.mismatches[item],
                gaps=self.gaps[item],
                q_start=self.q_starts[item],
                q_end=self.q_ends[item],
                t_start=self.t_starts[item],
                t_end=self.t_ends[item]
            )
        return PairwiseAlignments(
            scores=self.scores[item],
            matches=self.matches[item],
            mismatches=self.mismatches[item],
            gaps=self.gaps[item],
            q_starts=self.q_starts[item],
            q_ends=self.q_ends[item],
            t_starts=self.t_starts[item],
            t_ends=self.t_ends[item]
        )

    @classmethod
    def empty(cls) -> 'PairwiseAlignments':
        """Creates an empty PairwiseAlignments collection."""
        return cls(
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32)
        )

    @classmethod
    def concat(cls, batches: list['PairwiseAlignments']) -> 'PairwiseAlignments':
        """Concatenates multiple PairwiseAlignments collections into one."""
        if not batches:
            return cls.empty()
        return cls(
            scores=np.concatenate([b.scores for b in batches]),
            matches=np.concatenate([b.matches for b in batches]),
            mismatches=np.concatenate([b.mismatches for b in batches]),
            gaps=np.concatenate([b.gaps for b in batches]),
            q_starts=np.concatenate([b.q_starts for b in batches]),
            q_ends=np.concatenate([b.q_ends for b in batches]),
            t_starts=np.concatenate([b.t_starts for b in batches]),
            t_ends=np.concatenate([b.t_ends for b in batches]),
        )

    @property
    def pidents(self) -> npt.NDArray[np.float64]:
        """Calculates the percent identity for all alignments in the batch in a vectorized manner.

        Returns:
            npt.NDArray[np.float64]: A 1D array of percent identity values.
        """
        total = self.matches + self.mismatches + self.gaps
        return np.divide(
            self.matches * 100.0,
            total,
            out=np.zeros(len(self), dtype=np.float64),
            where=total > 0
        )


@dataclass(frozen=True, slots=True)
class PairwiseAligner:
    """A high-performance, batched pairwise sequence aligner.

    This class uses a parallelized, banded Smith-Waterman-Gotoh algorithm to align pairs of sequences.
    The core alignment logic is implemented in a Numba-jitted kernel for C-like speed. It is designed
    to operate on `Sequences` objects, aligning corresponding query and target sequences from two batches.

    The alignment can be guided by seeds (e.g., from a `RandstrobeIndex`) to constrain the dynamic
    programming matrix to a narrow band around a specific diagonal, drastically reducing computational
    complexity for closely related sequences.

    Attributes:
        gap_open (int): The penalty for opening a new gap. Defaults to 11.
        gap_extend (int): The penalty for extending an existing gap. Defaults to 1.
        k (int): The bandwidth parameter. The alignment is constrained to a band of `2*k+1` diagonals
            around the main diagonal (or a seed-defined diagonal). Defaults to 20.
    """
    gap_open: int = 11
    gap_extend: int = 1
    k: int = 20

    def __call__(self, queries: Sequences, targets: Sequences, seeds: 'Seeds | None' = None) -> PairwiseAlignments:
        """Performs pairwise alignment of query and target sequences.

        This method takes two `Sequences` objects of the same length and aligns each query sequence
        against its corresponding target sequence. The alignment uses the BLOSUM62 substitution matrix
        and the configured gap penalties.

        Args:
            queries (Sequences): A collection of query sequences.
            targets (Sequences): A collection of target sequences. Must have the same length as `queries`.
            seeds (Seeds, optional): Optional seeds to guide the alignment. If provided, the alignment
                for each sequence pair is centered on the diagonal specified by the corresponding seed.
                Defaults to None.

        Returns:
            PairwiseAlignments: A collection containing the scores and statistics for each alignment.

        Raises:
            ValueError: If the query and target batches have different numbers of sequences.
        """
        if len(queries.offsets) != len(targets.offsets):
            raise ValueError("Query and target batches must have the same number of sequences.")

        n = len(queries.offsets)
        if n == 0:
            return PairwiseAlignments.empty()

        # Handle the unified routing
        if seeds is not None:
            is_seeded = True
            offsets_arr = seeds.offsets
        else:
            is_seeded = False
            offsets_arr = np.zeros(n, dtype=np.int32)

        out_scores = np.empty(n, dtype=np.int32)
        out_matches = np.empty(n, dtype=np.int32)
        out_mismatches = np.empty(n, dtype=np.int32)
        out_gaps = np.empty(n, dtype=np.int32)
        out_q_starts = np.empty(n, dtype=np.int32)
        out_q_ends = np.empty(n, dtype=np.int32)
        out_t_starts = np.empty(n, dtype=np.int32)
        out_t_ends = np.empty(n, dtype=np.int32)
        
        _batched_banded_gotoh(
            queries.seqs, queries.offsets, queries.lengths,
            targets.seqs, targets.offsets, targets.lengths,
            _blosum62_matrix(), self.k, self.gap_open, self.gap_extend,
            is_seeded, offsets_arr,
            out_scores, out_matches, out_mismatches, out_gaps,
            out_q_starts, out_q_ends, out_t_starts, out_t_ends
        )
        
        return PairwiseAlignments(
            out_scores, out_matches, out_mismatches, out_gaps,
            out_q_starts, out_q_ends, out_t_starts, out_t_ends
        )

    def align_seeds(self, queries: Sequences, targets: Sequences, seeds: 'Seeds') -> PairwiseAlignments:
        """Convenience method to extract and align specific sequence pairs mapped by seeds.
        
        Args:
            queries (Sequences): The full collection of query sequences.
            targets (Sequences): The full collection of target sequences.
            seeds (Seeds): A collection of seeds mapping specific queries to targets.
            
        Returns:
            PairwiseAlignments: The resulting alignments, parallel to the provided Seeds.
        """
        paired_queries, paired_targets = seeds.extract_sequences(queries, targets)
        return self(paired_queries, paired_targets, seeds)


@dataclass(frozen=True, slots=True)
class Seeds:
    """A high-performance SoA container for alignment seeds.

    This class stores information about potential alignment regions (seeds) between query and target
    sequences. It is typically generated by a `RandstrobeIndex` search. The data is stored in a
    Structure-of-Arrays (SoA) layout for efficient processing.

    Attributes:
        query_indices (npt.NDArray[np.uint32]): A 1D array of indices specifying which query sequence
            each seed belongs to.
        target_indices (npt.NDArray[np.uint32]): A 1D array of indices specifying which target sequence
            each seed corresponds to.
        scores (npt.NDArray[np.uint32]): A 1D array of scores for each seed, typically representing
            the absolute count of shared randstrobes (the intersection size) between the sequences.
        offsets (npt.NDArray[np.int32]): A 1D array of diagonal offsets. The offset is calculated as
            `(query_pos - target_pos)` and is used to center the banded alignment.
    """
    query_indices: npt.NDArray[np.uint32]
    target_indices: npt.NDArray[np.uint32]
    scores: npt.NDArray[np.uint32]
    offsets: npt.NDArray[np.int32]
    
    def __len__(self) -> int:
        """Returns the number of seeds in the batch."""
        return len(self.query_indices)

    @classmethod
    def empty(cls) -> 'Seeds':
        """Creates an empty Seeds collection.
        
        Returns:
            Seeds: A new collection with zero-length arrays.
        """
        return cls(np.empty(0, dtype=np.uint32), np.empty(0, dtype=np.uint32),
                   np.empty(0, dtype=np.uint32), np.empty(0, dtype=np.int32))

    def filter(self, mask: npt.NDArray[np.bool_]) -> 'Seeds':
        """Returns a new Seeds collection containing only records where the mask is True.
        
        Args:
            mask (npt.NDArray[np.bool_]): A 1D boolean array used to filter the seeds.
            
        Returns:
            Seeds: A new, filtered collection.
        """
        return Seeds(
            self.query_indices[mask],
            self.target_indices[mask],
            self.scores[mask],
            self.offsets[mask]
        )

    def top_hits(self, min_score: int = 1) -> 'Seeds':
        """Reduces the batch to only the highest-scoring seed for each query.
        
        Args:
            min_score (int, optional): The minimum score for a seed to be returned.
                Defaults to 1, which filters out queries with zero score.
                
        Returns:
            Seeds: A new collection containing at most one seed per query.
        """
        if len(self) == 0:
            return self
            
        # np.lexsort prioritizes the last array. 
        # Primary key: query_indices (ascending). Secondary: scores (descending via ~)
        order = np.lexsort((~self.scores, self.query_indices))
        
        # Find the first occurrence of each query (which corresponds to its highest score)
        _, unique_idx = np.unique(self.query_indices[order], return_index=True)
        
        # Map back to the original array indices and sort to preserve the original query order
        best_idx = order[unique_idx]
        best_idx.sort()
        
        best_batch = Seeds(
            self.query_indices[best_idx],
            self.target_indices[best_idx],
            self.scores[best_idx],
            self.offsets[best_idx]
        )

        if min_score > 0:
            return best_batch.filter(best_batch.scores >= min_score)
            
        return best_batch

    def extract_sequences(self, queries: Sequences, targets: Sequences) -> tuple[Sequences, Sequences]:
        """Extracts parallel batches of query and target sequences mapped by this seed batch.

        Args:
            queries (Sequences): The source collection of query sequences.
            targets (Sequences): The source collection of target sequences.

        Returns:
            tuple[Sequences, Sequences]: Paired sequences ready for alignment.
        """
        return queries[self.query_indices], targets[self.target_indices]


@dataclass(frozen=True, slots=True)
class RandstrobeIndex:
    """A specialized index for fast, homology-based sequence search using randstrobes.

    This class implements an index based on the randstrobe method, which generates sparse,
    sequence-specific minimizers. Randstrobes are robust to local insertions and deletions,
    making them highly effective for finding homologous sequences.

    The index stores hashes of randstrobes and their corresponding positions. The `search` and
    `top_hits` methods perform a fast lookup of shared randstrobes between a query index and this
    (the target) index to identify candidate alignment regions (seeds).

    The core logic for building the index and performing searches is implemented in highly optimized,
    parallelized Numba kernels.

    Attributes:
        hashes (npt.NDArray[np.uint64]): A 1D array of sorted 64-bit randstrobe hash values.
        seq_indices (npt.NDArray[np.uint32]): A 1D array mapping each hash back to the index of the
            sequence it originated from in the source `Sequences` collection.
        pos1 (npt.NDArray[np.uint32]): The start position of the first k-mer (strobe 1) of the randstrobe.
        pos2 (npt.NDArray[np.uint32]): The start position of the second k-mer (strobe 2) of the randstrobe.
        n_seqs (int): The total number of sequences that were indexed.
        k (int): The k-mer size used for strobes.
        s (int): The s-mer size used for the syncmer selection within k-mers.
        w_min (int): The minimum window size for selecting the second strobe.
        w_max (int): The maximum window size for selecting the second strobe.
    """
    hashes:      npt.NDArray[np.uint64]
    seq_indices: npt.NDArray[np.uint32]
    pos1:        npt.NDArray[np.uint32]
    pos2:        npt.NDArray[np.uint32]
    n_seqs:      int = 0
    is_sorted:   bool = False
    k:           int = 10
    s:           int = 5
    w_min:       int = 1
    w_max:       int = 5

    def __len__(self) -> int:
        """Returns the total number of randstrobes in the index."""
        return len(self.hashes)
    
    @classmethod
    def empty(cls) -> 'RandstrobeIndex':
        """Creates an empty RandstrobeIndex."""
        return cls(
            np.empty(0, dtype=np.uint64),
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.uint32),
            0, False
        )

    @classmethod
    def build(cls, batch: Sequences, k: int = 10, s: int = 5, w_min: int = 1, w_max: int = 5, sort_by_hash: bool = False) -> 'RandstrobeIndex':
        """Builds a RandstrobeIndex from a `Sequences` collection.

        The build process involves two passes, both executed by parallel Numba kernels:

        1.  **Counting Pass:** Counts the number of randstrobes that will be generated from each
            sequence to pre-allocate the exact amount of memory needed.
        2.  **Population Pass:** Generates the randstrobes for each sequence and populates the
            pre-allocated NumPy arrays.

        Finally, if `sort_by_hash=True`, the arrays are sorted by hash value to enable fast binary searching during lookups.

        Args:
            batch (Sequences): The collection of sequences to index.
            k (int, optional): The k-mer size for strobes. Defaults to 10.
            s (int, optional): The s-mer (syncmer) size for selecting k-mers. Defaults to 5.
            w_min (int, optional): The minimum window offset for selecting the second strobe. Defaults to 1.
            w_max (int, optional): The maximum window offset for selecting the second strobe. Defaults to 5.
            sort_by_hash (bool, optional): Whether to sort the index by hash value. Required for target indices. Defaults to False.

        Returns:
            RandstrobeIndex: The constructed and sorted index.

        Raises:
            ValueError: If `s` is not strictly less than `k`.
        """

        if s >= k:
            raise ValueError("Sub-k-mer size (s) must be strictly less than k-mer size (k).")

        if len(batch) == 0:
            return cls.empty()

        # Pass 1: Count randstrobes per sequence to allocate exact memory footprint
        lut = _mmseqs12_lut()
        counts = _count_randstrobes_kernel(batch.seqs, batch.offsets, batch.lengths, lut, k, s, w_min)
        if (total_randstrobes := np.sum(counts)) == 0:
            return cls.empty()

        # Generate exact write offsets via cumulative sum
        out_offsets = np.empty(len(counts), dtype=np.uint32)
        current_offset = 0
        for i in range(len(counts)):
            out_offsets[i] = current_offset
            current_offset += counts[i]

        # Pre-allocate output arrays
        out_hashes = np.empty(total_randstrobes, dtype=np.uint64)
        out_seq_indices = np.empty(total_randstrobes, dtype=np.uint32)
        out_pos1 = np.empty(total_randstrobes, dtype=np.uint32)
        out_pos2 = np.empty(total_randstrobes, dtype=np.uint32)

        # Pass 2: Populate the randstrobes in parallel
        _populate_randstrobes_kernel(
            batch.seqs, batch.offsets, batch.lengths, lut, k, s, w_min, w_max, out_offsets,
            out_hashes, out_seq_indices, out_pos1, out_pos2
        )

        if sort_by_hash:
            sort_idx = np.argsort(out_hashes)
            out_hashes = out_hashes[sort_idx]
            out_seq_indices = out_seq_indices[sort_idx]
            out_pos1 = out_pos1[sort_idx]
            out_pos2 = out_pos2[sort_idx]

        return cls(
            out_hashes,
            out_seq_indices,
            out_pos1,
            out_pos2,
            len(batch), sort_by_hash, k, s, w_min, w_max
        )

    def _prep_queries(self, queries: 'RandstrobeIndex | Sequences') -> 'RandstrobeIndex':
        if not self.is_sorted:
            raise ValueError("Target index must be sorted by hash for binary search. Build it with sort_by_hash=True.")

        if isinstance(queries, Sequences):
            queries = RandstrobeIndex.build(queries, k=self.k, s=self.s, w_min=self.w_min, w_max=self.w_max)

        if queries.is_sorted:
            raise ValueError("Query index must NOT be sorted by hash. Build it with sort_by_hash=False.")

        return queries
    
    def top_hits(self, queries: 'RandstrobeIndex | Sequences', min_score: int = 1) -> Seeds:
        """Finds the single best-matching target sequence for each query sequence.

        This method iterates through each query sequence in the `queries` index, finds all shared
        randstrobe hashes with this (the target) index, and tallies the counts for each target sequence.
        The target with the highest number of shared hashes is reported as the top hit. By default,
        only hits with a score greater than zero are returned.

        Args:
            queries (RandstrobeIndex | Sequences): An index of query sequences, or `Sequences` to be
                indexed on the fly.
            min_score (int, optional): The minimum score for a seed to be returned. Defaults to 1, which
                filters out queries with no matching randstrobes. Set to 0 to include all queries.

        Returns:
            Seeds: A collection of seeds. If `min_score > 0`, the collection may contain fewer records
                than the number of input queries.
        """

        if len(queries) == 0 or len(self) == 0:
            return Seeds.empty()

        queries = self._prep_queries(queries)
        q_offsets = _compute_query_offsets(queries.seq_indices, queries.n_seqs)

        seeds = Seeds(*_intersect_top_hit_kernel(
            queries.hashes, queries.pos1, q_offsets, queries.n_seqs,
            self.hashes, self.seq_indices, self.pos1, self.n_seqs
        ))

        if min_score > 0:
            return seeds.filter(seeds.scores >= min_score)

        return seeds

    def search(self, queries: 'RandstrobeIndex | Sequences', min_score: int = 1, max_hits_per_query: int = 100) -> Seeds:
        """Finds all target sequences that meet a minimum score threshold for each query.

        Similar to `top_hits`, this method finds shared randstrobe hashes but reports all targets
        that have at least `min_score` shared hashes with a query.

        Args:
            queries (RandstrobeIndex): An index of query sequences.
            min_score (int, optional): The minimum number of shared hashes required to report a hit.
                Defaults to 1.
            max_hits_per_query (int, optional): The maximum number of hits to report for any single
                query. Defaults to 100.

        Returns:
            Seeds: A collection containing all seeds that passed the score threshold. The number of
                seeds can be greater than the number of queries.
        """
        if len(queries) == 0 or len(self) == 0:
            return Seeds.empty()

        queries = self._prep_queries(queries)
        q_offsets = _compute_query_offsets(queries.seq_indices, queries.n_seqs)

        # Pre-allocate 2D arrays for lock-free parallel writing
        out_t = np.zeros((queries.n_seqs, max_hits_per_query), dtype=np.uint32)
        out_s = np.zeros((queries.n_seqs, max_hits_per_query), dtype=np.uint32)
        out_o = np.zeros((queries.n_seqs, max_hits_per_query), dtype=np.int32)
        hit_counts = np.zeros(queries.n_seqs, dtype=np.uint32)

        _intersect_all_hits_kernel(
            queries.hashes, queries.pos1, q_offsets, queries.n_seqs,
            self.hashes, self.seq_indices, self.pos1, self.n_seqs,
            min_score, max_hits_per_query,
            out_t, out_s, out_o, hit_counts
        )

        # Generate repeated query indices based on how many hits each query got
        # Create a boolean mask to extract only the valid entries from the 2D arrays
        mask = np.arange(max_hits_per_query) < hit_counts[:, None]

        return Seeds(np.repeat(np.arange(queries.n_seqs, dtype=np.uint32), hit_counts),
                         out_t[mask], out_s[mask], out_o[mask])


# Functions ------------------------------------------------------------------------------------------------------------
@cache
def _blosum62_matrix(fill_value: int = -128) -> npt.NDArray[np.int8]:
    """Generates a 256x256-byte matrix for BLOSUM62 score lookup"""
    blosum62 = np.array([
        [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1, -1, -4],
        [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2, 0, -1, -4],
        [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 4, -3, 0, -1, -4],
        [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, -3, 1, -1, -4],
        [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4],
        [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, -2, 4, -1, -4],
        [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, -3, 4, -1, -4],
        [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -4],
        [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, -3, 0, -1, -4],
        [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, 3, -3, -1, -4],
        [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, 3, -3, -1, -4],
        [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, -3, 1, -1, -4],
        [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, 2, -1, -1, -4],
        [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, 0, -3, -1, -4],
        [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4],
        [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, -2, 0, -1, -4],
        [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1, -1, -4],
        [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -2, -2, -1, -4],
        [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -1, -2, -1, -4],
        [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, 2, -2, -1, -4],
        [-2, -1, 4, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, -3, 0, -1, -4],
        [-1, -2, -3, -3, -1, -2, -3, -4, -3, 3, 3, -3, 2, 0, -3, -2, -1, -2, -1, 2, -3, 3, -3, -1, -4],
        [-1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -2, -2, -2, 0, -3, 4, -1, -4],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4],
        [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]
    ], dtype=np.int8)
    alphabet = list(b"ARNDCQEGHILKMFPSTWYVBJZX*")
    score_matrix = np.full((256, 256), fill_value, dtype=np.int8)
    for a, i in enumerate(alphabet):
        for b, j in enumerate(alphabet):
            score_matrix[i, j] = blosum62[a, b]

    # Lock the array to make it strictly read-only and safe to share across instances
    score_matrix.flags.writeable = False
    return score_matrix


@cache
def _mmseqs12_lut(fill_value: int = 12) -> npt.NDArray[np.uint8]:
    """Generates a 256-byte array for O(1) MMSEQS12 alphabet mapping."""
    mapping = {
        b'A': 0, b'S': 0, b'T': 0,
        b'L': 1, b'M': 1,
        b'I': 2, b'V': 2,
        b'K': 3, b'R': 3,
        b'E': 4, b'Q': 4,
        b'N': 5, b'D': 5,
        b'F': 6, b'Y': 6,
        b'C': 7,
        b'G': 8,
        b'H': 9,
        b'P': 10,
        b'W': 11
    }
    # Map unknowns (e.g. X) to 12
    lut = np.full(256, fill_value, dtype=np.uint8)
    for source_byte, target_int in mapping.items():
        lut[source_byte[0]] = target_int
        lut[source_byte.lower()[0]] = target_int
    lut.flags.writeable = False
    return lut


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(parallel=True, cache=True, nogil=True)
def _batched_banded_gotoh(
        q_data: npt.NDArray[np.uint8], q_offsets: npt.NDArray[np.uint32], q_lengths: npt.NDArray[np.uint32],
        t_data: npt.NDArray[np.uint8], t_offsets: npt.NDArray[np.uint32], t_lengths: npt.NDArray[np.uint32],
        matrix: npt.NDArray[np.int8], k: int, gap_open: int, gap_extend: int,
        is_seeded: bool, diagonal_offsets: npt.NDArray[np.int32],
        out_scores: npt.NDArray[np.int32], out_matches: npt.NDArray[np.int32],
        out_mismatches: npt.NDArray[np.int32], out_gaps: npt.NDArray[np.int32],
        out_q_starts: npt.NDArray[np.int32], out_q_ends: npt.NDArray[np.int32],
        out_t_starts: npt.NDArray[np.int32], out_t_ends: npt.NDArray[np.int32]
):
    for idx in prange(len(q_offsets)):
        seq1 = q_data[q_offsets[idx]:q_offsets[idx] + q_lengths[idx]]
        seq2 = t_data[t_offsets[idx]:t_offsets[idx] + t_lengths[idx]]

        len1, len2 = len(seq1), len(seq2)
        rows, cols = len1 + 1, len2 + 1

        if is_seeded:  # Bandwidth is just the expected indels; offset shifts the center
            k_local = k
            offset = diagonal_offsets[idx]
        else:  # Classic Gotoh: band must absorb the sequence length difference; centered on 0
            k_local = max(k, abs(len1 - len2) + 1)
            offset = 0

        # TODO: Look into allocating the dynamic programming matrices as (rows, 2*k_local + 3) rather than (rows, cols),
        #  and addressing the columns using a sliding offset j_mapped = j - start_j.

        M = np.empty((rows, cols), dtype=np.int32)
        I = np.empty((rows, cols), dtype=np.int32)
        D = np.empty((rows, cols), dtype=np.int32)
        tb_M = np.empty((rows, cols), dtype=np.uint8)
        tb_D = np.empty((rows, cols), dtype=np.uint8)
        tb_I = np.empty((rows, cols), dtype=np.uint8)

        # O(N * k_local) band-only initialization
        for i in range(rows):
            j_center = i - offset
            start_j = max(0, j_center - k_local - 1)
            end_j = min(cols, j_center + k_local + 2)

            if start_j >= cols or end_j <= 0:
                continue

            for j in range(start_j, end_j):
                M[i, j] = 0
                I[i, j] = _INF
                D[i, j] = _INF
                tb_M[i, j] = 3  # 3 denotes the end of local alignment traceback

        max_score = 0
        max_i = 0
        max_j = 0

        # DP fill with local alignment logic (Smith-Waterman banded)
        for i in range(1, rows):
            j_center = i - offset
            start_j = max(1, j_center - k_local)
            end_j = min(cols, j_center + k_local + 1)

            if start_j >= cols or end_j <= 1:
                continue

            for j in range(start_j, end_j):
                d_open = M[i - 1, j] - gap_open - gap_extend
                d_ext = D[i - 1, j] - gap_extend
                if d_open >= d_ext:
                    D[i, j], tb_D[i, j] = d_open, 0
                else:
                    D[i, j], tb_D[i, j] = d_ext, 1

                i_open = M[i, j - 1] - gap_open - gap_extend
                i_ext = I[i, j - 1] - gap_extend
                if i_open >= i_ext:
                    I[i, j], tb_I[i, j] = i_open, 0
                else:
                    I[i, j], tb_I[i, j] = i_ext, 2

                m_diag = M[i - 1, j - 1] + matrix[seq1[i - 1], seq2[j - 1]]
                best, tb = m_diag, 0
                if D[i, j] > best:
                    best, tb = D[i, j], 1
                if I[i, j] > best:
                    best, tb = I[i, j], 2

                # Local alignment condition: score resets to 0 if negative
                if best <= 0:
                    M[i, j] = 0
                    tb_M[i, j] = 3
                else:
                    M[i, j] = best
                    tb_M[i, j] = tb
                    if best > max_score:
                        max_score = best
                        max_i = i
                        max_j = j

        # Traceback to count matches, mismatches, and gaps
        i, j = max_i, max_j
        matches = mismatches = gaps = state = 0

        # Store the end coordinates before traceback
        end_i, end_j = i, j

        while i > 0 and j > 0:
            if state == 0:
                tb = tb_M[i, j]
                if tb == 3:  # Reached the end of the local alignment
                    break
                elif tb == 0:  # Came from diagonal.
                    if seq1[i - 1] == seq2[j - 1]:
                        matches += 1
                    else:
                        mismatches += 1
                    i -= 1
                    j -= 1
                else:  # Came from D or I.
                    state = tb
            elif state == 1:  # Current state: in D. Execute an upward move.
                gaps += 1
                i -= 1
                if tb_D[i + 1, j] == 0:  # Gap open. We must have come from M.
                    state = 0
            elif state == 2:  # Current state: in I. Execute a leftward move.
                gaps += 1
                j -= 1
                if tb_I[i, j + 1] == 0:  # Gap open. We must have come from M.
                    state = 0

        # After loop, i and j are the start coordinates
        start_i, start_j = i, j

        out_scores[idx] = max_score
        out_matches[idx] = matches
        out_mismatches[idx] = mismatches
        out_gaps[idx] = gaps
        out_q_starts[idx] = start_i
        out_q_ends[idx] = end_i
        out_t_starts[idx] = start_j
        out_t_ends[idx] = end_j


@njit(cache=True, nogil=True, inline='always')
def _splitmix64(x: np.uint64) -> np.uint64:
    """Fast, invertible 64-bit integer mixing function for randomizing hashes."""
    x = np.uint64(x + np.uint64(0x9E3779B97F4A7C15))
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


@njit(parallel=True, cache=True, nogil=True)
def _count_randstrobes_kernel(
        seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.uint32], lengths: npt.NDArray[np.uint32],
        lut: npt.NDArray[np.uint8], k: int, s: int, w_min: int
) -> npt.NDArray[np.uint32]:
    n_seqs = len(offsets)
    counts = np.zeros(n_seqs, dtype=np.uint32)

    for idx in prange(n_seqs):
        seq_len = lengths[idx]
        if seq_len < k:
            continue

        seq_start = offsets[idx]
        n_syncmers = 0

        # Scan for random open syncmers
        for i in range(seq_len - k + 1):
            min_s_hash = ~np.uint64(0)
            min_s_idx = -1

            # Check all s-mers within this k-mer window
            for j in range(k - s + 1):
                s_val = np.uint64(0)
                for c in range(s):
                    char = seqs[seq_start + i + j + c]
                    s_val = s_val * np.uint64(12) + np.uint64(lut[char])

                h = _splitmix64(s_val)
                if h < min_s_hash:
                    min_s_hash = h
                    min_s_idx = j

            # Open Syncmer condition: minimum s-mer is at the first or last position
            if min_s_idx == 0 or min_s_idx == (k - s):
                n_syncmers += 1

        # A sequence yields exactly (n_syncmers - w_min) randstrobes (if positive)
        if n_syncmers > w_min:
            counts[idx] = n_syncmers - w_min

    return counts


@njit(parallel=True, cache=True, nogil=True)
def _populate_randstrobes_kernel(
        seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.uint32], lengths: npt.NDArray[np.uint32],
        lut: npt.NDArray[np.uint8], k: int, s: int, w_min: int, w_max: int, out_offsets: npt.NDArray[np.uint32],
        out_hashes: npt.NDArray[np.uint64], out_seq_indices: npt.NDArray[np.uint32],
        out_pos1: npt.NDArray[np.uint32], out_pos2: npt.NDArray[np.uint32]
):
    n_seqs = len(offsets)

    for idx in prange(n_seqs):
        seq_len = lengths[idx]
        if seq_len < k:
            continue

        seq_start = offsets[idx]
        max_syncmers = seq_len - k + 1

        # Thread-local allocations to store localized syncmer data safely within prange
        syncmer_pos = np.empty(max_syncmers, dtype=np.uint32)
        syncmer_hash = np.empty(max_syncmers, dtype=np.uint64)
        n_syncmers = 0

        # 1. Identify and cache the syncmers
        for i in range(seq_len - k + 1):
            min_s_hash = ~np.uint64(0)
            min_s_idx = -1

            for j in range(k - s + 1):
                s_val = np.uint64(0)
                for c in range(s):
                    char = seqs[seq_start + i + j + c]
                    s_val = s_val * np.uint64(12) + np.uint64(lut[char])

                h = _splitmix64(s_val)
                if h < min_s_hash:
                    min_s_hash = h
                    min_s_idx = j

            if min_s_idx == 0 or min_s_idx == (k - s):
                syncmer_pos[n_syncmers] = i

                # Hash the full k-mer window to use for randstrobe linkage
                k_val = np.uint64(0)
                for c in range(k):
                    char = seqs[seq_start + i + c]
                    k_val = k_val * np.uint64(12) + np.uint64(lut[char])

                syncmer_hash[n_syncmers] = _splitmix64(k_val)
                n_syncmers += 1

        # 2. Build Order-2 Randstrobes from the cached syncmers
        out_idx = out_offsets[idx]

        for i in range(n_syncmers - w_min):
            h1 = syncmer_hash[i]
            min_randstrobe_hash = ~np.uint64(0)
            best_j = -1

            end_j = i + w_max + 1
            if end_j > n_syncmers:
                end_j = n_syncmers

            for j in range(i + w_min, end_j):
                h2 = syncmer_hash[j]

                # Combine hashes asymmetrically (XOR of hashed values)
                combined = _splitmix64(h1 ^ _splitmix64(h2))
                if combined < min_randstrobe_hash:
                    min_randstrobe_hash = combined
                    best_j = j

            if best_j != -1:
                out_hashes[out_idx] = min_randstrobe_hash
                out_seq_indices[out_idx] = idx
                out_pos1[out_idx] = syncmer_pos[i]
                out_pos2[out_idx] = syncmer_pos[best_j]
                out_idx += 1


@njit(cache=True, nogil=True)
def _compute_query_offsets(q_seq_indices: npt.NDArray[np.uint32], num_queries: int) -> npt.NDArray[np.uint32]:
    """Finds the boundary indices for each query sequence in the flat hash array."""
    q_offsets = np.zeros(num_queries + 1, dtype=np.uint32)
    if len(q_seq_indices) == 0:
        return q_offsets

    current_query = 0
    for i in range(len(q_seq_indices)):
        q = q_seq_indices[i]
        while current_query < q:
            current_query += 1
            q_offsets[current_query] = i
            
    while current_query < num_queries:
        current_query += 1
        q_offsets[current_query] = len(q_seq_indices)
        
    return q_offsets


@njit(inline='always')
def _tally_single_query(
        start: int, end: int,
        q_hashes: npt.NDArray[np.uint64], q_pos: npt.NDArray[np.uint32],
        t_hashes: npt.NDArray[np.uint64], t_seq_idx: npt.NDArray[np.uint32], t_pos: npt.NDArray[np.uint32],
        tally: npt.NDArray[np.uint32], anchors: npt.NDArray[np.int32]
):
    for i in range(start, end):
        h = q_hashes[i]
        q_p = q_pos[i]

        # Find the first occurrence of the hash
        curr = np.searchsorted(t_hashes, h, side='left')

        # Walk forward to handle all collisions (handles out-of-bounds seamlessly)
        while curr < len(t_hashes) and t_hashes[curr] == h:
            t_id = t_seq_idx[curr]
            tally[t_id] += 1
            if tally[t_id] == 1:
                    anchors[t_id] = np.int32(q_p) - np.int32(t_pos[curr])
            curr += 1


@njit(parallel=True, cache=True, nogil=True)
def _intersect_top_hit_kernel(
        q_hashes, q_pos, q_offsets, num_queries,
        t_hashes, t_seq_idx, t_pos, num_targets
):
    best_t = np.zeros(num_queries, dtype=np.uint32)
    max_s = np.zeros(num_queries, dtype=np.uint32)
    offsets = np.zeros(num_queries, dtype=np.int32)

    for q_idx in prange(num_queries):
        start, end = q_offsets[q_idx], q_offsets[q_idx + 1]
        if start == end: continue

        # Thread-local allocations (Fits cleanly in L1 Cache)
        tally = np.zeros(num_targets, dtype=np.uint32)
        anchors = np.zeros(num_targets, dtype=np.int32)

        # The inline directive injects the binary search code directly here
        _tally_single_query(start, end, q_hashes, q_pos, t_hashes, t_seq_idx, t_pos, tally, anchors)

        # Argmax Search
        best_score, best_target = 0, 0
        for t in range(num_targets):
            if tally[t] > best_score:
                best_score = tally[t]
                best_target = t

        best_t[q_idx] = best_target
        max_s[q_idx] = best_score
        offsets[q_idx] = anchors[best_target]

    # Return dense arrays, returning q_idx explicitly for easy zipping
    q_idx_arr = np.arange(num_queries, dtype=np.uint32)
    return q_idx_arr, best_t, max_s, offsets


@njit(parallel=True, cache=True, nogil=True)
def _intersect_all_hits_kernel(
        q_hashes, q_pos, q_offsets, num_queries,
        t_hashes, t_seq_idx, t_pos, num_targets,
        min_score, max_hits_per_query,
        out_t, out_s, out_o, hit_counts
):
    for q_idx in prange(num_queries):
        start, end = q_offsets[q_idx], q_offsets[q_idx + 1]
        if start == end: continue

        tally = np.zeros(num_targets, dtype=np.uint32)
        anchors = np.zeros(num_targets, dtype=np.int32)

        # The identical inline code is injected here
        _tally_single_query(start, end, q_hashes, q_pos, t_hashes, t_seq_idx, t_pos, tally, anchors)

        # Matrix Population (Lock-Free)
        h_count = 0
        for t in range(num_targets):
            if tally[t] >= min_score:
                if h_count < max_hits_per_query:
                    out_t[q_idx, h_count] = t
                    out_s[q_idx, h_count] = tally[t]
                    out_o[q_idx, h_count] = anchors[t]
                    h_count += 1

        hit_counts[q_idx] = h_count
