from dataclasses import dataclass
from typing import Iterable, Self, NamedTuple
from enum import IntEnum

import numpy as np
import numpy.typing as npt
from numba import njit

from kaptive.core.constants import Strand
from kaptive.core.interval import IntervalBatch


# Classes --------------------------------------------------------------------------------------------------------------
class CigarOp(IntEnum):
    """BAM CIGAR operation encodings."""
    M = 0  # Alignment match (can be a sequence match or mismatch)
    I = 1  # Insertion to the reference
    D = 2  # Deletion from the reference
    N = 3  # Skipped region from the reference
    S = 4  # Soft clipping (clipped sequences present in SEQ)
    H = 5  # Hard clipping (clipped sequences NOT present in SEQ)
    P = 6  # Padding (silent deletion from padded reference)
    EQ = 7  # Sequence match
    X = 8  # Sequence mismatch
    B = 9  # Backwards compatibility

    @property
    def char(self) -> str:
        """Returns the string character for the operation."""
        return "MIDNSHP=XB"[self.value]


@dataclass(frozen=True, slots=True)
class CigarBatch:
    """Batched, flat-memory container for BAM-encoded CIGAR operations."""
    data: npt.NDArray[np.uint32]
    offsets: npt.NDArray[np.uint32]
    lengths: npt.NDArray[np.uint32]

    def __len__(self) -> int:
        return len(self.offsets)

    def __getitem__(self, item) -> npt.NDArray[np.uint32] | 'CigarBatch':
        """Access CIGAR arrays by index, slice, or boolean mask."""
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            s, l = self.offsets[item], self.lengths[item]
            return self.data[s:s + l]

        if isinstance(item, slice):
            indices = np.arange(len(self))[item]
        else:
            item_arr = np.asarray(item)
            indices = np.nonzero(item_arr)[0] if item_arr.dtype.kind == 'b' else item_arr

        if len(indices) == 0:
            return self.empty()

        new_lengths = self.lengths[indices]
        new_offsets = np.zeros(len(new_lengths), dtype=np.uint32)
        if len(new_lengths) > 1:
            np.cumsum(new_lengths[:-1], out=new_offsets[1:])

        extracted = np.concatenate([self.data[self.offsets[i]:self.offsets[i] + self.lengths[i]] for i in indices])
        return CigarBatch(extracted, new_offsets, new_lengths)

    @classmethod
    def empty(cls) -> 'CigarBatch':
        return cls(np.empty(0, dtype=np.uint32), np.empty(0, dtype=np.uint32), np.empty(0, dtype=np.uint32))

    @classmethod
    def concat(cls, batches: Iterable['CigarBatch']) -> 'CigarBatch':
        batches = list(batches)
        if not batches:
            return cls.empty()
        lengths = np.concatenate([b.lengths for b in batches])
        offsets = np.zeros(len(lengths), dtype=np.uint32)
        if len(lengths) > 1:
            np.cumsum(lengths[:-1], out=offsets[1:])
        return cls(np.concatenate([b.data for b in batches]), offsets, lengths)

    def swap_sides(self) -> 'CigarBatch':
        """Returns a new CigarBatch with Insertions and Deletions flipped."""
        return CigarBatch(_swap_cigar_kernel(self.data), self.offsets, self.lengths)

    @classmethod
    def from_lists(cls, cigar_lists: list[npt.NDArray[np.uint32]]) -> 'CigarBatch':
        if not cigar_lists:
            return cls.empty()
        lengths = np.array([len(c) for c in cigar_lists], dtype=np.uint32)
        offsets = np.zeros(len(lengths), dtype=np.uint32)
        if len(lengths) > 1:
            np.cumsum(lengths[:-1], out=offsets[1:])
        return cls(np.concatenate(cigar_lists), offsets, lengths)


class AlignmentRecord(NamedTuple):
    """
    A lightweight, self-aware view of a single alignment record.

    Attributes:
        idx: Original index within the source AlignmentBatch.
        q_name: Query sequence name.
        q_length: Total length of the query sequence.
        q_start: Start position on the query.
        q_end: End position on the query.
        t_name: Target (reference) sequence name.
        t_length: Total length of the target sequence.
        t_start: Start position on the target.
        t_end: End position on the target.
        strand: Alignment orientation (Strand.FORWARD or Strand.REVERSE).
        length: Total alignment block length.
        match: Number of matching bases.
        mismatch: Number of mismatches (NM tag).
        quality: Mapping quality (MAPQ).
        cigar: CIGAR array
    """
    idx: int
    q_name: str
    q_length: int
    q_start: int
    q_end: int
    t_name: str
    t_length: int
    t_start: int
    t_end: int
    strand: Strand
    length: int
    match: int
    mismatch: int
    quality: int
    cigar: npt.NDArray[np.uint32]


@dataclass(frozen=True, slots=True)
class AlignmentBatch:
    """
    A high-performance batch of alignments stored in NumPy arrays.

    Uses a Structure-of-Arrays (SoA) layout to enable fast vectorized operations
    and efficient filtering of large alignment sets.
    """
    q_names: npt.NDArray[np.object_]
    q_lengths: npt.NDArray[np.int32]
    q_starts: npt.NDArray[np.int32]
    q_ends: npt.NDArray[np.int32]
    t_names: npt.NDArray[np.object_]
    t_lengths: npt.NDArray[np.int32]
    t_starts: npt.NDArray[np.int32]
    t_ends: npt.NDArray[np.int32]
    strands: npt.NDArray[np.int8]
    lengths: npt.NDArray[np.int32]
    matches: npt.NDArray[np.int32]
    mismatches: npt.NDArray[np.int32]
    qualities: npt.NDArray[np.int32]
    cigars: CigarBatch

    def __len__(self) -> int:
        return len(self.q_starts)

    @property
    def scores(self) -> np.ndarray:
        """Vectorized alignment score (matches - mismatches)."""
        return self.matches - self.mismatches

    @classmethod
    def from_mappy(cls, q_name: str, q_length: int, alignments: Iterable['mappy.Alignment']) -> Self:
        """
        Create an AlignmentBatch from mappy Alignment objects.

        Args:
            q_name: Name of the query sequence.
            q_length: Length of the query sequence.
            alignments: Iterable of mappy.Alignment objects.

        Returns:
            AlignmentBatch: A new batch containing the alignments.
        """
        qn, ql, qs, qe, tn, tl, ts, te, st, bl, ml, nm, mq = [], [], [], [], [], [], [], [], [], [], [], [], []
        cigar_lists = []

        for h in alignments:
            qn.append(q_name)
            ql.append(q_length)
            qs.append(h.q_st)
            qe.append(h.q_en)
            tn.append(h.ctg)
            tl.append(h.ctg_len)
            ts.append(h.r_st)
            te.append(h.r_en)
            st.append(h.strand)
            bl.append(h.blen)
            ml.append(h.mlen)
            nm.append(h.NM)
            mq.append(h.mapq)
            cigar_lists.append(np.array([(l << 4) | op for l, op in h.cigar], dtype=np.uint32))

        if not qn:
            raise ValueError("Cannot initialize AlignmentBatch with empty alignments")

        return cls(
            q_names=np.array(qn, dtype=object),
            q_lengths=np.array(ql, dtype=np.int32),
            q_starts=np.array(qs, dtype=np.int32),
            q_ends=np.array(qe, dtype=np.int32),
            t_names=np.array(tn, dtype=object),
            t_lengths=np.array(tl, dtype=np.int32),
            t_starts=np.array(ts, dtype=np.int32),
            t_ends=np.array(te, dtype=np.int32),
            strands=np.array(st, dtype=np.int8),
            lengths=np.array(bl, dtype=np.int32),
            matches=np.array(ml, dtype=np.int32),
            mismatches=np.array(nm, dtype=np.int32),
            qualities=np.array(mq, dtype=np.int32),
            cigars=CigarBatch.from_lists(cigar_lists)
        )

    @classmethod
    def concat(cls, batches: Iterable['AlignmentBatch']) -> Self:
        """Concatenate multiple AlignmentBatch objects into one."""
        batches = list(batches)
        if not batches:
            raise ValueError("Cannot concatenate an empty iterable of batches")

        kwargs = {}
        for field_name in cls.__dataclass_fields__:
            if field_name == 'cigars':
                kwargs[field_name] = CigarBatch.concat([b.cigars for b in batches])
                continue
            first_val = getattr(batches[0], field_name)
            if isinstance(first_val, np.ndarray):
                kwargs[field_name] = np.concatenate([getattr(b, field_name) for b in batches])
            else:
                if any(getattr(b, field_name) != first_val for b in batches):
                    raise ValueError(f"Cannot concatenate batches with mismatched '{field_name}' values")
                kwargs[field_name] = first_val

        return cls(**kwargs)

    def __getitem__(self, item) -> 'AlignmentRecord | AlignmentBatch':
        """Access alignments by index, slice, or boolean mask."""
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return AlignmentRecord(
                idx=item, q_name=self.q_names[item], q_length=self.q_lengths[item],
                q_start=self.q_starts[item], q_end=self.q_ends[item], t_name=self.t_names[item],
                t_length=self.t_lengths[item], t_start=self.t_starts[item], t_end=self.t_ends[item],
                strand=self.strands[item], length=self.lengths[item], match=self.matches[item],
                mismatch=self.mismatches[item], quality=self.qualities[item], cigar=self.cigars[item]
            )
            
        return AlignmentBatch(
            q_names=self.q_names[item], q_lengths=self.q_lengths[item], q_starts=self.q_starts[item],
            q_ends=self.q_ends[item], t_names=self.t_names[item], t_lengths=self.t_lengths[item],
            t_starts=self.t_starts[item], t_ends=self.t_ends[item], strands=self.strands[item],
            lengths=self.lengths[item], matches=self.matches[item], mismatches=self.mismatches[item],
            qualities=self.qualities[item], cigars=self.cigars[item]
        )

    def filter(self, mask: np.ndarray) -> 'AlignmentBatch':
        """Return a new batch containing only elements where mask is True. (Alias for self[mask])"""
        return self[mask]

    def filter_out(self, mask: np.ndarray) -> 'AlignmentBatch':
        """Returns a new batch excluding the masked items."""
        return self.filter(~mask)

    def cull_overlaps(self, max_overlap_fraction: float = 0.1, group_by: np.ndarray | None = None,
                      priority_mask: np.ndarray | None = None) -> 'AlignmentBatch':
        """
        Greedily culls overlapping alignments on the target using a fast Numba kernel.
        """
        if (n := len(self)) < 2:
            return self
        
        # Map string t_names to integer IDs for C-level kernel processing
        _, t_name_ints = np.unique(self.t_names, return_inverse=True)
        
        if priority_mask is not None:
            sort_scores = self.scores.astype(np.float64)
            sort_scores[priority_mask] += 1e9
            order = np.argsort(sort_scores, kind='stable')[::-1].astype(np.int32)
        else:
            order = np.argsort(self.scores, kind='stable')[::-1].astype(np.int32)
        
        if group_by is None:
            group_by = np.zeros(n, dtype=np.int32)
        
        kept_mask = _cull_overlaps_kernel(
            order, t_name_ints, self.t_starts, self.t_ends, group_by, max_overlap_fraction, n
        )
        return self.filter(kept_mask)

    def shift_query(self, offset: int) -> 'AlignmentBatch':
        """
        Returns a new AlignmentBatch with query coordinates shifted by the given offset.
        Used when mapping chunks of a sequence in parallel.
        """
        if offset == 0 or len(self) == 0:
            return self
        
        return AlignmentBatch(
            q_names=self.q_names, q_lengths=self.q_lengths, 
            q_starts=self.q_starts + offset, q_ends=self.q_ends + offset, 
            t_names=self.t_names, t_lengths=self.t_lengths,
            t_starts=self.t_starts, t_ends=self.t_ends, strands=self.strands,
            lengths=self.lengths, matches=self.matches, mismatches=self.mismatches,
            qualities=self.qualities, cigars=self.cigars
        )

    def swap_sides(self) -> 'AlignmentBatch':
        """
        Swaps query and target roles in the alignment records.

        This is used when query sequences (contigs) are being treated as targets
        and vice-versa, which is common in reciprocal mapping.
        """
        return AlignmentBatch(
            q_names=self.t_names,
            q_lengths=self.t_lengths,
            q_starts=self.t_starts,
            q_ends=self.t_ends,
            t_names=self.q_names,
            t_lengths=self.q_lengths,
            t_starts=self.q_starts,
            t_ends=self.q_ends,
            strands=self.strands,
            lengths=self.lengths,
            matches=self.matches,
            mismatches=self.mismatches,
            qualities=self.qualities,
            cigars=self.cigars.swap_sides()
        )

    def split(self, by_query: bool = False) -> Iterable[tuple[str, 'AlignmentBatch']]:
        """Splits a batch into separate batches by target or query name."""
        key_array = self.q_names if by_query else self.t_names
        for key in np.unique(key_array):
            yield key, self.filter(key_array == key)

    def to_intervals(self, by_query: bool = False) -> IntervalBatch:
        """
        Converts the batch into a high-performance IntervalBatch.

        Args:
            by_query: If True, use query coordinates. Otherwise, use target coordinates.
        """
        starts = self.q_starts if by_query else self.t_starts
        ends = self.q_ends if by_query else self.t_ends

        return IntervalBatch(
            starts=starts,
            ends=ends,
            strands=self.strands,
            # CRITICAL: Ensures we can map relational queries back to this alignment record!
            original_indices=np.arange(len(self), dtype=np.int32)
        )

    def get_record(self, idx: int) -> AlignmentRecord:
        """Retrieve a single AlignmentRecord by its batch index. (Alias for self[idx])"""
        return self[idx]

    @property
    def is_partial_mask(self) -> np.ndarray:
        """Vectorized boolean mask of alignments covering < 90% of the query."""
        return self.lengths < (self.q_lengths * 0.9)

    @classmethod
    def from_records(cls, records: Iterable[AlignmentRecord]) -> 'AlignmentBatch':
        """
        Create an AlignmentBatch from an iterable of AlignmentRecord objects.
        
        Uses a fast list comprehension and zip transposition to efficiently 
        reconstruct the Structure of Arrays (SoA) layout.
        """
        data = [
            (r.q_name, r.q_length, r.q_start, r.q_end, r.t_name, r.t_length,
             r.t_start, r.t_end, r.strand, r.length, r.match, r.mismatch,
             r.quality)
            for r in records

        ]
        if not data:
            raise ValueError("Cannot initialize AlignmentBatch with empty records")
        qn, ql, qs, qe, tn, tl, ts, te, st, bl, ml, nm, mq = zip(*data)
        return cls(
            q_names=np.array(qn, dtype=object),
            q_lengths=np.array(ql, dtype=np.int32),
            q_starts=np.array(qs, dtype=np.int32),
            q_ends=np.array(qe, dtype=np.int32),
            t_names=np.array(tn, dtype=object),
            t_lengths=np.array(tl, dtype=np.int32),
            t_starts=np.array(ts, dtype=np.int32),
            t_ends=np.array(te, dtype=np.int32),
            strands=np.array(st, dtype=np.int8),
            lengths=np.array(bl, dtype=np.int32),
            matches=np.array(ml, dtype=np.int32),
            mismatches=np.array(nm, dtype=np.int32),
            qualities=np.array(mq, dtype=np.int32),
            cigars=CigarBatch.from_lists([r.cigar for r in records])
        )


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(cache=True, nogil=True)
def _cull_overlaps_kernel(order: npt.NDArray[np.int32], t_names: npt.NDArray[np.int32], starts: npt.NDArray[np.int32],
                          ends: npt.NDArray[np.int32], groups: npt.NDArray[np.int32], max_overlap_fraction: float, n: int):
    """Numba-accelerated greedy overlap culling."""
    kept_mask = np.zeros(n, dtype=np.bool_)
    
    for i in range(n):
        idx = order[i]
        t = t_names[idx]
        s = starts[idx]
        e = ends[idx]
        g = groups[idx]
        length = e - s
        
        if length <= 0:
            continue
            
        overlap_found = False
        # Check against previously kept intervals in O(N^2) (Very fast in pure C)
        for j in range(i):
            prev_idx = order[j]
            if not kept_mask[prev_idx] or t_names[prev_idx] != t or groups[prev_idx] != g:
                continue
                
            ks, ke = starts[prev_idx], ends[prev_idx]
            overlap = min(e, ke) - max(s, ks)
            if overlap > 0 and (overlap / length) > max_overlap_fraction:
                overlap_found = True
                break
                
        if not overlap_found:
            kept_mask[idx] = True
            
    return kept_mask


@njit(cache=True, nogil=True)
def parse_cigar_string(cigar_str: str) -> npt.NDArray[np.uint32]:
    """Fast Numba parser converting CIGAR string to BAM-encoded uint32 array."""
    ops = 0
    for char in cigar_str:
        if char == 'M' or char == 'I' or char == 'D' or char == 'N' or char == 'S' or char == 'H' or char == 'P' or char == '=' or char == 'X' or char == 'B':
            ops += 1

    out = np.empty(ops, dtype=np.uint32)
    idx = 0
    current_len = 0

    for char in cigar_str:
        if '0' <= char <= '9':
            current_len = current_len * 10 + (ord(char) - ord('0'))
        else:
            op = 0
            if char == 'M':
                op = 0
            elif char == 'I':
                op = 1
            elif char == 'D':
                op = 2
            elif char == 'N':
                op = 3
            elif char == 'S':
                op = 4
            elif char == 'H':
                op = 5
            elif char == 'P':
                op = 6
            elif char == '=':
                op = 7
            elif char == 'X':
                op = 8
            elif char == 'B':
                op = 9

            out[idx] = (np.uint32(current_len) << 4) | np.uint32(op)
            idx += 1
            current_len = 0

    return out


@njit(cache=True, nogil=True)
def _swap_cigar_kernel(cigar_data: npt.NDArray[np.uint32]) -> npt.NDArray[np.uint32]:
    """Swaps Insertions (1) with Deletions (2) at the bit-level."""
    out = np.empty_like(cigar_data)
    for i in range(len(cigar_data)):
        val = cigar_data[i]
        op = val & 0xF
        if op == 1:
            out[i] = (val & ~np.uint32(0xF)) | np.uint32(2)
        elif op == 2:
            out[i] = (val & ~np.uint32(0xF)) | np.uint32(1)
        else:
            out[i] = val
    return out