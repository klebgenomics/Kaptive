r"""Pairwise sequence alignment algorithms and containers.

This module provides high-performance pairwise sequence alignment capabilities
using a parallelized, banded Smith-Waterman-Gotoh dynamic programming algorithm
with affine gap penalties. Alignment results are stored in memory-efficient Structure-of-Arrays
(SoA) containers ([`PairwiseAlignments`][kaptive.core.pairwise.PairwiseAlignments]) or
scalar records ([`PairwiseAlignment`][kaptive.core.pairwise.PairwiseAlignment]).

Key Classes:
    - [`PairwiseAlignment`][kaptive.core.pairwise.PairwiseAlignment]: Immutable container for a single result.
    - [`PairwiseAlignments`][kaptive.core.pairwise.PairwiseAlignments]: Batched SoA container for results.
    - [`PairwiseAligner`][kaptive.core.pairwise.PairwiseAligner]: Banded Smith-Waterman-Gotoh alignment engine.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from functools import cache
from typing import Any, Self

import numpy as np
import numpy.typing as npt
from numba import njit, prange

from kaptive.core.collections import BatchedContainer
from kaptive.core.kmers import Seeds
from kaptive.core.seq import Sequences

# Constants ------------------------------------------------------------------------------------------------------------
_INF = np.int32(-1_000_000_000)  # Safe negative infinity for int32 DP arithmetic


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class PairwiseAlignment:
    r"""A lightweight, immutable container for the results of a single pairwise sequence alignment.

    This class holds summary statistics and coordinates for an alignment between a single query
    and target sequence. It is typically produced by indexing into a
    [`PairwiseAlignments`][kaptive.core.pairwise.PairwiseAlignments] collection.

    Attributes:
        score (int): Final alignment score calculated using BLOSUM62 and gap penalties.
        matches (int): Total number of matching bases.
        mismatches (int): Total number of mismatched bases.
        gaps (int): Total number of gap characters (insertions or deletions).
        q_start (int): 0-based start coordinate on query sequence (inclusive).
        q_end (int): 0-based end coordinate on query sequence (exclusive).
        t_start (int): 0-based start coordinate on target sequence (inclusive).
        t_end (int): 0-based end coordinate on target sequence (exclusive).
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
        r"""Calculate the percent identity of the aligned region.

        Percent identity is defined as `(matches / (matches + mismatches + gaps)) * 100.0`.

        Returns:
            float: Percent identity value between 0.0 and 100.0, or 0.0 if alignment length is zero.
        """
        total = self.matches + self.mismatches + self.gaps
        return (self.matches / total) * 100.0 if total > 0 else 0.0


@dataclass(frozen=True, slots=True)
class PairwiseAlignments(BatchedContainer["PairwiseAlignment", "PairwiseAlignments"]):
    r"""A high-performance SoA container for the results of multiple pairwise alignments.

    This class stores alignment statistics in a Structure-of-Arrays (SoA) layout using 1D NumPy arrays.

    Attributes:
        scores (npt.NDArray[np.int32]): 1D array of alignment scores.
        matches (npt.NDArray[np.int32]): 1D array of match counts.
        mismatches (npt.NDArray[np.int32]): 1D array of mismatch counts.
        gaps (npt.NDArray[np.int32]): 1D array of gap counts.
        q_starts (npt.NDArray[np.int32]): 1D array of query start coordinates.
        q_ends (npt.NDArray[np.int32]): 1D array of query end coordinates.
        t_starts (npt.NDArray[np.int32]): 1D array of target start coordinates.
        t_ends (npt.NDArray[np.int32]): 1D array of target end coordinates.
    """

    scores: npt.NDArray[np.int32]
    matches: npt.NDArray[np.int32]
    mismatches: npt.NDArray[np.int32]
    gaps: npt.NDArray[np.int32]
    q_starts: npt.NDArray[np.int32]
    q_ends: npt.NDArray[np.int32]
    t_starts: npt.NDArray[np.int32]
    t_ends: npt.NDArray[np.int32]

    def __len__(self) -> int:
        r"""Return the number of alignments in the batch.

        Returns:
            int: Total count of alignments.
        """
        return len(self.scores)

    def to_dict(self) -> dict[str, npt.NDArray[np.int32]]:
        r"""Convert the alignment batch to a dictionary of NumPy arrays for serialization.

        Returns:
            dict[str, npt.NDArray[np.int32]]: Dictionary mapping attribute names to arrays.
        """
        return {
            "scores": self.scores,
            "matches": self.matches,
            "mismatches": self.mismatches,
            "gaps": self.gaps,
            "q_starts": self.q_starts,
            "q_ends": self.q_ends,
            "t_starts": self.t_starts,
            "t_ends": self.t_ends,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> PairwiseAlignments:
        r"""Deserialize a PairwiseAlignments batch from a dictionary of array-like data.

        Args:
            d (dict[str, Any]): Dictionary containing alignment field arrays.

        Returns:
            PairwiseAlignments: Deserialized pairwise alignments container.
        """
        return cls(
            np.array(d["scores"], dtype=np.int32),
            np.array(d["matches"], dtype=np.int32),
            np.array(d["mismatches"], dtype=np.int32),
            np.array(d["gaps"], dtype=np.int32),
            np.array(d["q_starts"], dtype=np.int32),
            np.array(d["q_ends"], dtype=np.int32),
            np.array(d["t_starts"], dtype=np.int32),
            np.array(d["t_ends"], dtype=np.int32),
        )

    def __getitem__(self, item: Any) -> PairwiseAlignment | PairwiseAlignments:
        r"""Access alignment results by index, slice, or boolean array mask.

        Args:
            item (Any): Integer index, slice, or boolean/integer array.

        Returns:
            PairwiseAlignment | PairwiseAlignments: Scalar record or sliced batch collection.

        Raises:
            IndexError: If integer index is out of range.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return PairwiseAlignment(
                score=int(self.scores[item]),
                matches=int(self.matches[item]),
                mismatches=int(self.mismatches[item]),
                gaps=int(self.gaps[item]),
                q_start=int(self.q_starts[item]),
                q_end=int(self.q_ends[item]),
                t_start=int(self.t_starts[item]),
                t_end=int(self.t_ends[item]),
            )
        return PairwiseAlignments(
            scores=self.scores[item],
            matches=self.matches[item],
            mismatches=self.mismatches[item],
            gaps=self.gaps[item],
            q_starts=self.q_starts[item],
            q_ends=self.q_ends[item],
            t_starts=self.t_starts[item],
            t_ends=self.t_ends[item],
        )

    @classmethod
    def empty(cls) -> PairwiseAlignments:
        r"""Create an empty PairwiseAlignments collection with zero-length int32 arrays.

        Returns:
            PairwiseAlignments: Empty alignment collection.
        """
        return cls(
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenate multiple PairwiseAlignments collections into a single batch.

        Args:
            batches (Iterable[PairwiseAlignments]): Iterable of alignment collections.

        Returns:
            PairwiseAlignments: Single concatenated alignment collection.
        """
        batches_list = list(batches)
        if not batches_list:
            return cls.empty()  # type: ignore
        return cls(
            scores=np.concatenate([b.scores for b in batches_list]),
            matches=np.concatenate([b.matches for b in batches_list]),
            mismatches=np.concatenate([b.mismatches for b in batches_list]),
            gaps=np.concatenate([b.gaps for b in batches_list]),
            q_starts=np.concatenate([b.q_starts for b in batches_list]),
            q_ends=np.concatenate([b.q_ends for b in batches_list]),
            t_starts=np.concatenate([b.t_starts for b in batches_list]),
            t_ends=np.concatenate([b.t_ends for b in batches_list]),
        )

    @property
    def pidents(self) -> npt.NDArray[np.float64]:
        r"""Calculate percent identity for all alignments in the batch in a vectorized manner.

        Returns:
            npt.NDArray[np.float64]: 1D array of percent identity values.
        """
        total = self.matches + self.mismatches + self.gaps
        return np.divide(self.matches * 100.0, total, out=np.zeros(len(self), dtype=np.float64), where=total > 0)


@dataclass(frozen=True, slots=True)
class PairwiseAligner:
    r"""A high-performance, batched pairwise sequence aligner.

    Uses a parallelized, banded Smith-Waterman-Gotoh algorithm to align pairs of sequences.

    Attributes:
        gap_open (int): Penalty for opening a new gap. Defaults to 11.
        gap_extend (int): Penalty for extending an existing gap. Defaults to 1.
        k (int): Bandwidth parameter (`2*k+1` diagonals). Defaults to 20.
    """

    gap_open: int = 11
    gap_extend: int = 1
    k: int = 20

    def __call__(self, queries: Sequences, targets: Sequences, seeds: Seeds | None = None) -> PairwiseAlignments:
        r"""Perform pairwise alignment of query and target sequences.

        Args:
            queries (Sequences): Collection of query sequences.
            targets (Sequences): Collection of target sequences (must match query batch size).
            seeds (Seeds | None): Optional alignment seeds guiding diagonal alignment.

        Returns:
            PairwiseAlignments: Alignment scores, statistics, and coordinates for each sequence pair.

        Raises:
            ValueError: If query and target batches have different numbers of sequences.
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
            queries.seqs,
            queries.offsets,
            queries.lengths,
            targets.seqs,
            targets.offsets,
            targets.lengths,
            _blosum62_matrix(),
            self.k,
            self.gap_open,
            self.gap_extend,
            is_seeded,
            offsets_arr,
            out_scores,
            out_matches,
            out_mismatches,
            out_gaps,
            out_q_starts,
            out_q_ends,
            out_t_starts,
            out_t_ends,
        )

        return PairwiseAlignments(
            out_scores,
            out_matches,
            out_mismatches,
            out_gaps,
            out_q_starts,
            out_q_ends,
            out_t_starts,
            out_t_ends,
        )

    def align_seeds(self, queries: Sequences, targets: Sequences, seeds: Seeds) -> PairwiseAlignments:
        r"""Convenience method to extract and align specific sequence pairs mapped by seeds.

        Args:
            queries (Sequences): Full collection of query sequences.
            targets (Sequences): Full collection of target sequences.
            seeds (Seeds): Seed collection mapping specific query sequences to target sequences.

        Returns:
            PairwiseAlignments: Alignment results parallel to the provided seeds.
        """
        paired_queries, paired_targets = seeds.extract_sequences(queries, targets)
        return self(paired_queries, paired_targets, seeds)


# Functions ------------------------------------------------------------------------------------------------------------
@cache
def _blosum62_matrix(fill_value: int = -128) -> npt.NDArray[np.int8]:
    r"""Generate a 256x256-byte matrix for BLOSUM62 score lookup.

    Args:
        fill_value (int): Default score value for non-standard characters. Defaults to -128.

    Returns:
        npt.NDArray[np.int8]: Read-only 256x256 8-bit integer substitution matrix array.
    """
    blosum62 = np.array(
        [
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
            [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1],
        ],
        dtype=np.int8,
    )
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
    q_data: npt.NDArray[np.uint8],
    q_offsets: npt.NDArray[np.int32],
    q_lengths: npt.NDArray[np.int32],
    t_data: npt.NDArray[np.uint8],
    t_offsets: npt.NDArray[np.int32],
    t_lengths: npt.NDArray[np.int32],
    matrix: npt.NDArray[np.int8],
    k: int,
    gap_open: int,
    gap_extend: int,
    is_seeded: bool,
    diagonal_offsets: npt.NDArray[np.int32],
    out_scores: npt.NDArray[np.int32],
    out_matches: npt.NDArray[np.int32],
    out_mismatches: npt.NDArray[np.int32],
    out_gaps: npt.NDArray[np.int32],
    out_q_starts: npt.NDArray[np.int32],
    out_q_ends: npt.NDArray[np.int32],
    out_t_starts: npt.NDArray[np.int32],
    out_t_ends: npt.NDArray[np.int32],
) -> None:
    r"""Numba parallel kernel for batched banded Smith-Waterman-Gotoh pairwise alignment.

    Args:
        q_data (npt.NDArray[np.uint8]): Concatenated query sequence bytes.
        q_offsets (npt.NDArray[np.int32]): Query sequence start offsets.
        q_lengths (npt.NDArray[np.int32]): Query sequence lengths.
        t_data (npt.NDArray[np.uint8]): Concatenated target sequence bytes.
        t_offsets (npt.NDArray[np.int32]): Target sequence start offsets.
        t_lengths (npt.NDArray[np.int32]): Target sequence lengths.
        matrix (npt.NDArray[np.int8]): 256x256 BLOSUM62 substitution score matrix.
        k (int): Alignment bandwidth parameter.
        gap_open (int): Gap opening penalty.
        gap_extend (int): Gap extension penalty.
        is_seeded (bool): Whether alignment is guided by diagonal seed offsets.
        diagonal_offsets (npt.NDArray[np.int32]): Diagonal offsets per sequence pair when seeded.
        out_scores (npt.NDArray[np.int32]): Output array for alignment scores.
        out_matches (npt.NDArray[np.int32]): Output array for match counts.
        out_mismatches (npt.NDArray[np.int32]): Output array for mismatch counts.
        out_gaps (npt.NDArray[np.int32]): Output array for gap counts.
        out_q_starts (npt.NDArray[np.int32]): Output array for query start coordinates.
        out_q_ends (npt.NDArray[np.int32]): Output array for query end coordinates.
        out_t_starts (npt.NDArray[np.int32]): Output array for target start coordinates.
        out_t_ends (npt.NDArray[np.int32]): Output array for target end coordinates.
    """
    for idx in prange(len(q_offsets)):  # type: ignore
        seq1 = q_data[q_offsets[idx] : q_offsets[idx] + q_lengths[idx]]
        seq2 = t_data[t_offsets[idx] : t_offsets[idx] + t_lengths[idx]]

        len1, len2 = len(seq1), len(seq2)
        rows, cols = len1 + 1, len2 + 1

        if is_seeded:  # Bandwidth is just the expected indels; offset shifts the center
            k_local = k
            offset = diagonal_offsets[idx]
        else:  # Classic Gotoh: band must absorb the sequence length difference; centered on 0
            k_local = max(k, abs(len1 - len2) + 1)
            offset = 0

        # Allocate dynamic programming matrices as (rows, band_width) to save huge amounts of memory
        band_width = 2 * k_local + 3
        M = np.empty((rows, band_width), dtype=np.int32)
        I_mat = np.empty((rows, band_width), dtype=np.int32)  # noqa: E741
        D = np.empty((rows, band_width), dtype=np.int32)
        tb_M = np.empty((rows, band_width), dtype=np.uint8)
        tb_D = np.empty((rows, band_width), dtype=np.uint8)
        tb_I = np.empty((rows, band_width), dtype=np.uint8)

        # O(N * k_local) band-only initialization
        for i in range(rows):
            j_center = i - offset
            start_j = max(0, j_center - k_local - 1)
            end_j = min(cols, j_center + k_local + 2)

            if start_j >= cols or end_j <= 0:
                continue

            for j in range(start_j, end_j):
                jm = j - start_j
                M[i, jm] = 0
                I_mat[i, jm] = _INF
                D[i, jm] = _INF
                tb_M[i, jm] = 3  # 3 denotes the end of local alignment traceback

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

            start_j_prev = max(0, i - 1 - offset - k_local - 1)
            start_j_curr = max(0, i - offset - k_local - 1)

            for j in range(start_j, end_j):
                jm_top = j - start_j_prev
                d_open = M[i - 1, jm_top] - gap_open - gap_extend
                d_ext = D[i - 1, jm_top] - gap_extend

                jm = j - start_j_curr
                if d_open >= d_ext:
                    D[i, jm], tb_D[i, jm] = d_open, 0
                else:
                    D[i, jm], tb_D[i, jm] = d_ext, 1

                jm_left = j - 1 - start_j_curr
                i_open = M[i, jm_left] - gap_open - gap_extend
                i_ext = I_mat[i, jm_left] - gap_extend
                if i_open >= i_ext:
                    I_mat[i, jm], tb_I[i, jm] = i_open, 0
                else:
                    I_mat[i, jm], tb_I[i, jm] = i_ext, 2

                jm_top_left = j - 1 - start_j_prev
                m_diag = M[i - 1, jm_top_left] + matrix[seq1[i - 1], seq2[j - 1]]

                best, tb = m_diag, 0
                if D[i, jm] > best:
                    best, tb = D[i, jm], 1
                if I_mat[i, jm] > best:
                    best, tb = I_mat[i, jm], 2

                # Local alignment condition: score resets to 0 if negative
                if best <= 0:
                    M[i, jm] = 0
                    tb_M[i, jm] = 3
                else:
                    M[i, jm] = best
                    tb_M[i, jm] = tb
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
            start_j_curr = max(0, i - offset - k_local - 1)
            jm = j - start_j_curr

            if state == 0:
                tb = tb_M[i, jm]
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
                tb = tb_D[i, jm]
                gaps += 1
                i -= 1
                if tb == 0:  # Gap open. We must have come from M.
                    state = 0
            elif state == 2:  # Current state: in I. Execute a leftward move.
                tb = tb_I[i, jm]
                gaps += 1
                j -= 1
                if tb == 0:  # Gap open. We must have come from M.
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
