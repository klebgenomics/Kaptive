from dataclasses import dataclass
from kaptive.core.seq import SeqBatch

import numpy as np
import numpy.typing as npt
from numba import njit, prange


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class PairwiseAlignment:
    score: int
    matches: int
    mismatches: int
    gaps: int

    @property
    def pident(self) -> float:
        """Returns the percent identity of the aligned region."""
        total = self.matches + self.mismatches + self.gaps
        return (self.matches / total) * 100.0 if total > 0 else 0.0


@dataclass(frozen=True, slots=True)
class PairwiseAlignmentBatch:
    scores:     npt.NDArray[np.int32]
    matches:    npt.NDArray[np.int32]
    mismatches: npt.NDArray[np.int32]
    gaps:       npt.NDArray[np.int32]

    def __len__(self) -> int:
        return len(self.scores)

    def __getitem__(self, item) -> 'PairwiseAlignment | PairwiseAlignmentBatch':
        """Access alignments by index, slice, or boolean mask."""
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return PairwiseAlignment(
                score=self.scores[item],
                matches=self.matches[item],
                mismatches=self.mismatches[item],
                gaps=self.gaps[item]
            )
        return PairwiseAlignmentBatch(
            scores=self.scores[item],
            matches=self.matches[item],
            mismatches=self.mismatches[item],
            gaps=self.gaps[item]
        )

    @classmethod
    def empty(cls) -> 'PairwiseAlignmentBatch':
        return cls(
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32)
        )

    @property
    def pidents(self) -> npt.NDArray[np.float64]:
        """Vectorized percent identity calculation."""
        total = self.matches + self.mismatches + self.gaps
        pidents = np.zeros(len(self), dtype=np.float64)
        valid = total > 0
        pidents[valid] = (self.matches[valid] / total[valid]) * 100.0
        return pidents


class PairwiseAligner:
    """
    Class for aligning pairs of amino-acid sequences using the banded gotoh algorithm.
    """
    __slots__ = ('_matrix', '_k', '_gap_open', '_gap_extend') # Use __slots__ for memory efficiency
    def __init__(self, gap_open: int = 11, gap_extend: int = 1, k: int = 20):
        self._gap_open = gap_open
        self._gap_extend = gap_extend
        self._k = k
        self._matrix = _blosum62()

    def __call__(self, queries: SeqBatch, targets: SeqBatch) -> PairwiseAlignmentBatch:
        if len(queries.offsets) != len(targets.offsets):
            raise ValueError("Query and target batches must have the same number of sequences.")

        n = len(queries.offsets)
        if n == 0:
            return PairwiseAlignmentBatch.empty()

        # Pre-allocate array for Numba to fill in parallel
        results_arr = np.empty((n, 4), dtype=np.int32)

        _batched_banded_gotoh(
            queries.seqs, queries.offsets, queries.lengths,
            targets.seqs, targets.offsets, targets.lengths,
            self._matrix, self._k, self._gap_open, self._gap_extend,
            results_arr
        )

        return PairwiseAlignmentBatch(
            scores=results_arr[:, 0],
            matches=results_arr[:, 1],
            mismatches=results_arr[:, 2],
            gaps=results_arr[:, 3]
        )


# Functions ------------------------------------------------------------------------------------------------------------
def _blosum62(empty_val: int = -128) -> npt.NDArray[np.int8]:
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
    score_matrix = np.full((256, 256), empty_val, dtype=np.int8)
    for a, i in enumerate(alphabet):
        for b, j in enumerate(alphabet):
            score_matrix[i, j] = blosum62[a, b]

    # Lock the array to make it strictly read-only and safe to share across instances
    score_matrix.flags.writeable = False
    return score_matrix


# Constants ------------------------------------------------------------------------------------------------------------
_INF = np.int32(-1_000_000_000)  # Safe negative infinity for int32 DP arithmetic


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(parallel=True, cache=True, nogil=True)
def _batched_banded_gotoh(
        q_data: npt.NDArray[np.uint8], q_offsets: npt.NDArray[np.uint32], q_lengths: npt.NDArray[np.uint32],
        t_data: npt.NDArray[np.uint8], t_offsets: npt.NDArray[np.uint32], t_lengths: npt.NDArray[np.uint32],
        matrix: npt.NDArray[np.int8], k, gap_open: int, gap_extend: int,
        results: npt.NDArray[np.int32]
):
    for idx in prange(len(q_offsets)):
        seq1 = q_data[q_offsets[idx]:q_offsets[idx] + q_lengths[idx]]
        seq2 = t_data[t_offsets[idx]:t_offsets[idx] + t_lengths[idx]]

        len1, len2 = len(seq1), len(seq2)
        rows, cols = len1 + 1, len2 + 1
        k_local = max(k, abs(len1 - len2) + 1)

        M = np.empty((rows, cols), dtype=np.int32)
        I = np.empty((rows, cols), dtype=np.int32)
        D = np.empty((rows, cols), dtype=np.int32)
        tb_M = np.empty((rows, cols), dtype=np.uint8)
        tb_D = np.empty((rows, cols), dtype=np.uint8)
        tb_I = np.empty((rows, cols), dtype=np.uint8)

        # O(N * k) band-only initialisation — this is the entire point of banding
        for i in range(rows):
            start_j = max(0, i - k_local - 1)
            end_j = min(cols, i + k_local + 2)
            for j in range(start_j, end_j):
                M[i, j] = _INF
                I[i, j] = _INF
                D[i, j] = _INF

        M[0, 0] = 0

        for i in range(1, min(rows, k_local + 1)):
            D[i, 0] = -gap_open - i * gap_extend  # fixed: was (i - 1)
            tb_D[i, 0] = 0 if i == 1 else 1
            M[i, 0] = D[i, 0]
            tb_M[i, 0] = 1

        for j in range(1, min(cols, k_local + 1)):
            I[0, j] = -gap_open - j * gap_extend  # fixed: was (j - 1)
            tb_I[0, j] = 0 if j == 1 else 2
            M[0, j] = I[0, j]
            tb_M[0, j] = 2

        # DP fill (unchanged — already correctly banded)
        for i in range(1, rows):
            start_j = max(1, i - k_local)
            end_j = min(cols, i + k_local + 1)
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
                if D[i, j] > best: best, tb = D[i, j], 1
                if I[i, j] > best: best, tb = I[i, j], 2
                M[i, j], tb_M[i, j] = best, tb

        # --- Traceback to count matches, mismatches, and gaps ---
        i, j = len1, len2
        matches = mismatches = gaps = 0
        # The traceback is a state machine. 'state' determines the current move.
        # 0: Process a diagonal move (from M).
        # 1: Process an upward move (from D).
        # 2: Process a leftward move (from I).
        # Initialize the state with the best move for the final cell.
        state = tb_M[i, j]
        while i > 0 or j > 0:
            if state == 0:  # Current state: in M. Decide next move.
                if i > 0 and j > 0:
                    # Check the pointer of the current M cell to see where it came from.
                    tb = tb_M[i, j]
                    if tb == 0:  # Came from diagonal.
                        # Execute the diagonal move.
                        if seq1[i - 1] == seq2[j - 1]: matches += 1
                        else: mismatches += 1
                        i -= 1; j -= 1
                        # The next state is still M, so state remains 0.
                    else:  # Came from D or I.
                        # Switch to that state for the next iteration.
                        state = tb
                elif i > 0: state = 1  # Edge case: must be deletions
                elif j > 0: state = 2  # Edge case: must be insertions
                else: break
            elif state == 1:  # Current state: in D. Execute an upward move.
                if i == 0: break
                gaps += 1; i -= 1
                # After moving, check if we opened or extended the gap.
                if tb_D[i + 1, j] == 0:  # Gap open. We must have come from M.
                    state = 0  # Switch back to M state for the next iteration.
                # else: Gap extension. Stay in D state. state remains 1.
            elif state == 2:  # Current state: in I. Execute a leftward move.
                if j == 0: break
                gaps += 1; j -= 1
                # After moving, check if we opened or extended the gap.
                if tb_I[i, j + 1] == 0:  # Gap open. We must have come from M.
                    state = 0  # Switch back to M state for the next iteration.
                # else: Gap extension. Stay in I state. state remains 2.

        score = M[len1, len2]
        # --- End Inlined Kernel ---

        results[idx, 0] = score
        results[idx, 1] = matches
        results[idx, 2] = mismatches
        results[idx, 3] = gaps
