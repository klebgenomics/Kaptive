"""
Module for performing pairwise protein alignments using a numba-powered banded Gotoh Smith-Waterman kernel.
"""
from dataclasses import dataclass
from functools import cache

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

    def __getitem__(self, item) -> 'PairwiseAlignment | PairwiseAlignments':
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

