"""
Genomic interval representation with strand and context, plus batched interval operations.
"""
from collections.abc import Iterable
from dataclasses import dataclass
from re import Match
from typing import Self, Union

import numpy as np
import numpy.typing as npt
from numba import njit

from kaptive.core.constants import Context, Strand


# Classes --------------------------------------------------------------------------------------------------------------
IntervalLike = Union[slice, int, Match, 'Interval']  # Type alias


@dataclass(frozen=True, slots=True)
class Interval:
    """
    A single genomic interval defined by start, end, and strand.
    Uses 0-based coordinate system (start inclusive, end exclusive).
    """
    start: int
    end: int
    strand: Strand = Strand.UNSTRANDED

    def __len__(self) -> int:
        return self.end - self.start

    def __contains__(self, item: IntervalLike) -> bool:
        """Check if an coordinate or another interval is fully contained within this one."""
        if isinstance(item, int): return self.start <= item < self.end
        item = Interval.from_item(item)
        return self.start <= item.start and self.end >= item.end

    def __add__(self, other: IntervalLike) -> 'Interval':
        """Returns the minimal bounding interval covering both self and other."""
        other = Interval.from_item(other)
        new_strand = self.strand if self.strand == other.strand else Strand.UNSTRANDED
        return Interval(min(self.start, other.start), max(self.end, other.end), new_strand)

    def __radd__(self, other: IntervalLike) -> 'Interval':
        return self.__add__(other)

    def shift(self, x: int, y: int | None = None) -> 'Interval':
        """
        Shift the interval by a fixed distance.

        Args:
            x: Distance to shift the start.
            y: Optional distance to shift the end. If None, same as x.
        """
        return Interval(self.start + x, self.end + (y if y is not None else x), self.strand)

    def expand(self, left: int, right: int, clip_length: int | None = None) -> 'Interval':
        """
        Expand (or shrink) the interval by specified amounts on the left and right.
        Safely clips the start coordinate to 0, and optionally clips the end coordinate.
        """
        new_start = max(0, self.start - left)
        new_end = self.end + right
        if clip_length is not None:
            new_end = min(new_end, clip_length)
            
        return Interval(new_start, new_end, self.strand)

    def reverse_complement(self, length: int | None = None) -> 'Interval':
        """
        Returns the interval coordinates on the opposite strand.

        Args:
            length: The length of the parent sequence (e.g. contig).
        """
        if length is None:
            length = self.end
        return Interval(length - self.end, length - self.start, self.strand * -1)

    def relate(self, other: IntervalLike) -> Context:
        """
        Calculate the spatial relationship between this interval and another.

        Returns a Context enum value (e.g. UPSTREAM, INSIDE, OVERLAPPING).
        """
        other = Interval.from_item(other)
        return Context(_core_relate(self.start, self.end, self.strand, other.start, other.end))

    @classmethod
    def from_match(cls, item: Match, strand: Strand = Strand.UNSTRANDED) -> 'Interval':
        """Create from a regex Match object."""
        return cls(item.start(), item.end(), strand)

    @classmethod
    def from_int(cls, item: int, strand: Strand = Strand.UNSTRANDED, length: int | None = None) -> 'Interval':
        """Create a 1bp interval from an integer."""
        if item < 0 and length is not None: item += length
        return cls(item, item + 1, strand)

    @classmethod
    def from_slice(cls, item: slice, strand: Strand = Strand.UNSTRANDED, length: int | None = None) -> 'Interval':
        """Create from a Python slice object."""
        start, stop, step = item.start, item.stop, item.step
        if start is None:
            start = 0
        if stop is None and length is not None:
            stop = length
        if stop is None:
            raise ValueError("Cannot create Interval from slice with None stop without 'length'")
        if step == -1:
            return cls(stop + 1, start + 1, strand)
        return cls(start, stop, strand)

    @classmethod
    def from_item(cls, item: IntervalLike, strand: Strand = Strand.UNSTRANDED, length: int | None = None) -> 'Interval':
        """Universal coercion from various objects to an Interval."""
        if isinstance(item, cls):
            return item
        if interval := getattr(item, 'interval', None) is not None:
            return interval
        if isinstance(item, Match):
            return cls.from_match(item, strand)
        if isinstance(item, int):
            return cls.from_int(item, strand, length)
        if isinstance(item, slice):
            return cls.from_slice(item, strand, length)
        raise TypeError(item)


@dataclass(frozen=True, slots=True)
class IntervalBatch:
    """
    High-performance batch of genomic intervals, powered by NumPy.

    Uses strict Structure-of-Arrays (SoA) layout with automatic dtype enforcement.
    Optimized for spatial queries (overlap, nearest-neighbor) and vectorized
    coordinate transformations.
    """
    starts: npt.NDArray[np.int32]
    ends: npt.NDArray[np.int32]
    strands: npt.NDArray[np.int8]
    original_indices: npt.NDArray[np.int32] | None = None

    def __post_init__(self):
        if self.original_indices is None:
            object.__setattr__(self, 'original_indices', np.arange(len(self.starts), dtype=np.int32))

    @classmethod
    def empty(cls) -> Self:
        """Create an empty IntervalBatch."""
        return cls(
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8),
            np.empty(0, dtype=np.int32)
        )

    @classmethod
    def from_intervals(cls, intervals: Iterable[Interval]) -> Self:
        """Create a batch from an iterable of Interval objects."""
        # OPTIMIZATION: Fast C-level list comprehension + zip extraction
        data = [(i.start, i.end, i.strand) for i in intervals]
        if not data:
            return cls.empty()

        s, e, st = zip(*data, strict=False)
        return cls(
            np.array(s, dtype=np.int32),
            np.array(e, dtype=np.int32),
            np.array(st, dtype=np.int8)
        )

    def max_len(self) -> int:
        """Returns the length of the longest interval in the batch."""
        return np.max(self.ends - self.starts) if len(self) > 0 else 0

    def __len__(self) -> int:
        return len(self.starts)

    def __getitem__(self, item):
        """
        Access intervals by index.

        If index is an integer, returns a single Interval object.
        If index is a slice or mask, returns a new IntervalBatch.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return Interval(self.starts[item], self.ends[item], self.strands[item])

        return IntervalBatch(
            self.starts[item],
            self.ends[item],
            self.strands[item],
            self.original_indices[item]
        )

    def copy(self):
        """Create a deep copy of the batch arrays."""
        return IntervalBatch(
            self.starts.copy(),
            self.ends.copy(),
            self.strands.copy(),
            self.original_indices.copy()
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:
        """Concatenate multiple IntervalBatch objects."""
        batches = list(batches)
        if not batches:
            raise ValueError("Cannot concatenate empty list of batches")
        return cls(
            np.concatenate([b.starts for b in batches]),
            np.concatenate([b.ends for b in batches]),
            np.concatenate([b.strands for b in batches]),
            np.concatenate([b.original_indices for b in batches])
        )

    def sort(self) -> 'IntervalBatch':
        """Returns a new cleanly sorted IntervalBatch (by start then end)."""
        if len(self) < 2 or _is_sorted_kernel(self.starts, self.ends):
            return self

        order = np.lexsort((self.ends, self.starts))
        return IntervalBatch(
            self.starts[order],
            self.ends[order],
            self.strands[order],
            self.original_indices[order]
        )

    def filter(self, mask: IntervalLike) -> 'IntervalBatch':
        """Return a new batch containing only intervals matched by the mask."""
        if isinstance(mask, (slice, int, np.integer)):
            if isinstance(mask, (int, np.integer)):
                mask = [mask]
            return self[mask]
        return self[np.asarray(mask)]

    @property
    def centers(self) -> npt.NDArray:
        """Vectorized calculation of interval midpoints."""
        return (self.starts + self.ends) / 2

    @property
    def lengths(self) -> npt.NDArray[np.int32]:
        """Vectorized calculation of interval lengths."""
        return self.ends - self.starts

    def query(self, start: int, end: int) -> npt.NDArray[np.int32]:
        """
        Find indices of intervals that overlap with the given range.

        Args:
            start: Query range start.
            end: Query range end.

        Returns:
            np.ndarray: Indices of overlapping intervals.
        """
        if len(self) == 0:
            return np.empty(0, dtype=self._INDEX_DTYPE)
        return _query_kernel(self.starts, self.ends, start, end, self.max_len())

    def merge(self, tolerance: int = 0) -> 'IntervalBatch':
        """
        Merge overlapping or adjacent intervals into single bounding boxes.

        Args:
            tolerance: Max gap between intervals to consider them adjacent.
        """
        if len(self) == 0:
            return self
        out = _merge_kernel(self.starts, self.ends, self.strands, tolerance)
        return type(self)(out[0], out[1], out[2])

    def expand(self, left: int | npt.NDArray, right: int | npt.NDArray, 
               clip_lengths: int | npt.NDArray | None = None) -> 'IntervalBatch':
        """
        Expand (or shrink) intervals by specified amounts on the left and right.
        Safely clips start coordinates to 0, and optionally clips end coordinates.
        """
        new_starts = np.maximum(0, self.starts - left)
        new_ends = self.ends + right
        
        if clip_lengths is not None:
            new_ends = np.minimum(new_ends, clip_lengths)
            
        return IntervalBatch(
            starts=new_starts.astype(self.starts.dtype),
            ends=new_ends.astype(self.ends.dtype),
            strands=self.strands.copy(),
            original_indices=self.original_indices.copy()
        )

    def project(self, shift: int, flip_length: int | None = None) -> 'IntervalBatch':
        """
        Apply a coordinate transformation to all intervals in the batch.

        Args:
            shift: Distance to add to all coordinates.
            flip_length: If provided, coordinates are mirrored within this length
                (e.g. for projecting onto the reverse strand of a contig).
        """
        if flip_length is not None:
            new_starts = (flip_length - self.ends).astype(self._START_DTYPE)
            new_ends = (flip_length - self.starts).astype(self._END_DTYPE)
            new_strands = (-self.strands).astype(self._STRAND_DTYPE)
        else:
            new_starts = self.starts.copy()
            new_ends = self.ends.copy()
            new_strands = self.strands.copy()

        new_starts += shift
        new_ends += shift

        return IntervalBatch(
            starts=new_starts,
            ends=new_ends,
            strands=new_strands,
            original_indices=self.original_indices.copy()
        )

    def overlaps_with(self, other: 'IntervalBatch') -> np.ndarray:
        """
        Vectorized many-to-many overlap check using 2D broadcasting.
        Returns a boolean mask of intervals in `self` that overlap ANY interval in `other`.
        """
        if len(self) == 0 or len(other) == 0:
            return np.zeros(len(self), dtype=bool)

        overlaps = (self.starts[:, None] < other.ends[None, :]) & \
                   (other.starts[None, :] < self.ends[:, None])
        return overlaps.any(axis=1)


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(cache=True, nogil=True)
def _query_kernel(starts: npt.NDArray[np.int32], ends: npt.NDArray[np.int32], q_start: int, q_end: int, max_len: int):
    """Numba-accelerated spatial overlap query."""
    limit = np.searchsorted(starts, q_end, side='left')
    count = 0
    min_start_check = q_start - max_len

    for i in range(limit - 1, -1, -1):
        if starts[i] < min_start_check:
            break
        if ends[i] > q_start:
            count += 1

    out = np.empty(count, dtype=np.int32)
    idx = 0
    for i in range(limit - 1, -1, -1):
        if starts[i] < min_start_check:
            break
        if ends[i] > q_start:
            out[idx] = i
            idx += 1

    return out[::-1]

@njit(cache=True, nogil=True)
def _merge_kernel(starts: npt.NDArray[np.int32], ends: npt.NDArray[np.int32], strands: npt.NDArray[np.int8],
                  tolerance: int) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.int8]]:
    """Numba-accelerated interval merging."""
    n = len(starts)
    if n == 0:
        return (np.empty(0, dtype=starts.dtype),
                np.empty(0, dtype=ends.dtype),
                np.empty(0, dtype=strands.dtype))

    temp_s = np.empty(n, dtype=starts.dtype)
    temp_e = np.empty(n, dtype=ends.dtype)
    temp_st = np.empty(n, dtype=strands.dtype)

    curr_s = starts[0]
    curr_e = ends[0]
    curr_st = strands[0]
    out_idx = 0

    for i in range(1, n):
        s = starts[i]
        e = ends[i]
        st = strands[i]

        if s <= curr_e + tolerance:
            curr_e = max(curr_e, e)
            if curr_st != st:
                curr_st = 0
        else:
            temp_s[out_idx] = curr_s
            temp_e[out_idx] = curr_e
            temp_st[out_idx] = curr_st
            out_idx += 1

            curr_s = s
            curr_e = e
            curr_st = st

    temp_s[out_idx] = curr_s
    temp_e[out_idx] = curr_e
    temp_st[out_idx] = curr_st
    out_idx += 1

    return temp_s[:out_idx], temp_e[:out_idx], temp_st[:out_idx]


@njit(cache=True, nogil=True)
def _is_sorted_kernel(starts: npt.NDArray[np.int32], ends: npt.NDArray[np.int32]) -> bool:
    """Numba-accelerated check for sortedness."""
    n = len(starts)
    for i in range(n - 1):
        if starts[i] > starts[i + 1]:
            return False
        if starts[i] == starts[i + 1] and ends[i] > ends[i + 1]:
            return False
    return True


@njit(cache=True, nogil=True, inline='always')
def _core_relate(s_a, e_a, st_a, s_b, e_b) -> int:
    """Core logic for determining spatial relationship between two intervals."""
    if e_b <= s_a:
        return 1 if st_a >= 0 else 2
    if s_b >= e_a:
        return 2 if st_a >= 0 else 1
    if s_b >= s_a and e_b <= e_a:
        return 3
    if s_b < s_a:
        if e_b > e_a:
            return 4
        return 5 if st_a >= 0 else 6
    return 6 if st_a >= 0 else 5
