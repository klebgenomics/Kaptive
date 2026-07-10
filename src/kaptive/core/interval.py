"""
Genomic interval representation with strand and context, plus batched interval operations.
"""
from collections.abc import Iterable
from dataclasses import dataclass
from enum import IntEnum
from re import Match
from typing import Union

import numpy as np
import numpy.typing as npt
from numba import njit


# Enums ----------------------------------------------------------------------------------------------------------------
class Strand(IntEnum):
    """
    Integer representation of genomic strand orientation.

    Supports conversion from common string formats (+, -, 1, -1).
    """
    FORWARD = 1
    REVERSE = -1
    UNSTRANDED = 0

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, bytes):
            value = value.decode('ascii')
        if isinstance(value, str):
            if value == '+' or value == '1' or value == '+1':
                return Strand.FORWARD
            elif value == '-' or value == '-1':
                return Strand.REVERSE
        return Strand.UNSTRANDED

    def __str__(self):
        if self == Strand.FORWARD:
            return '+'
        if self == Strand.REVERSE:
            return '-'
        return '.'


# Classes --------------------------------------------------------------------------------------------------------------
IntervalLike = Union[slice, int, Match, 'Interval']  # Type alias


@dataclass(frozen=True, slots=True)
class Interval:
    """A single genomic interval defined by start, end, and strand.

    This class represents a discrete segment on a biological sequence, such as a gene on a chromosome.
    It uses a 0-based, half-open coordinate system, where the `start` position is inclusive and the `end`
    position is exclusive. This is consistent with Python's slicing and common bioinformatics formats.

    The `strand` attribute indicates the orientation on the sequence, using the `Strand` enum
    (FORWARD: 1, REVERSE: -1, UNSTRANDED: 0).

    The class is implemented as a frozen dataclass with slots for high performance and immutability.
    It supports basic interval arithmetic (containment, union) and geometric transformations
    (shifting, expanding, reverse complementing).

    Attributes:
        start (int): The 0-based starting position of the interval (inclusive).
        end (int): The 0-based ending position of the interval (exclusive).
        strand (Strand): The orientation of the interval on its parent sequence.
    """
    start: int
    end: int
    strand: Strand = Strand.UNSTRANDED

    def __len__(self) -> int:
        """Calculates the length of the interval (end - start)."""
        return self.end - self.start

    def __contains__(self, item: IntervalLike) -> bool:
        """Checks if a coordinate or another interval is fully contained within this one.

        Args:
            item (IntervalLike): An integer coordinate or another interval-like object
                (Interval, slice, regex Match).

        Returns:
            bool: True if the item is completely inside the interval's bounds.
        """
        if isinstance(item, int): return self.start <= item < self.end
        item = Interval.from_item(item)
        return self.start <= item.start and self.end >= item.end

    def __add__(self, other: IntervalLike) -> 'Interval':
        """Computes the minimal bounding interval that covers both this interval and another.

        This is equivalent to finding the union of the two intervals' ranges. The strand of the
        resulting interval is preserved only if both input intervals have the same strand; otherwise,
        it becomes UNSTRANDED.

        Args:
            other (IntervalLike): The interval to merge with this one.

        Returns:
            Interval: A new Interval representing the combined span.
        """
        other = Interval.from_item(other)
        new_strand = self.strand if self.strand == other.strand else Strand.UNSTRANDED
        return Interval(min(self.start, other.start), max(self.end, other.end), new_strand)

    def __radd__(self, other: IntervalLike) -> 'Interval':
        """Reverse addition to support `other + self` syntax."""
        return self.__add__(other)

    def shift(self, x: int, y: int | None = None) -> 'Interval':
        """Shifts the interval by a fixed distance.

        This is useful for converting between coordinate systems (e.g., from local gene coordinates
        to global contig coordinates).

        Args:
            x (int): The distance to shift the start position.
            y (int, optional): The distance to shift the end position. If None, the end is shifted
                by the same amount as the start (`x`). Defaults to None.

        Returns:
            Interval: A new Interval with the shifted coordinates.
        """
        return Interval(self.start + x, self.end + (y if y is not None else x), self.strand)

    def expand(self, left: int, right: int, clip_length: int | None = None) -> 'Interval':
        """Expands or shrinks the interval by specified amounts on each side.

        This method is useful for creating flanking regions or padding. The start coordinate is
        automatically clipped at 0 to prevent negative values. The end coordinate can optionally
        be clipped to a maximum length.

        Args:
            left (int): The amount to expand the start (subtract from `start`). Can be negative to shrink.
            right (int): The amount to expand the end (add to `end`). Can be negative to shrink.
            clip_length (int, optional): The maximum length of the parent sequence. If provided,
                the new end coordinate will not exceed this value. Defaults to None.

        Returns:
            Interval: A new, expanded Interval.
        """
        new_start = max(0, self.start - left)
        new_end = self.end + right
        if clip_length is not None:
            new_end = min(new_end, clip_length)
            
        return Interval(new_start, new_end, self.strand)

    def reverse_complement(self, length: int | None = None) -> 'Interval':
        """Calculates the coordinates of the interval on the opposite strand.

        This operation mirrors the interval's position relative to a parent sequence of a given length
        and inverts its strand. For example, an interval `(10, 20)` on a sequence of length 100 would
        become `(80, 90)` on the reverse strand.

        Args:
            length (int, optional): The length of the parent sequence (e.g., a contig). If not provided,
                the interval's own end is used, which is suitable for intervals at the start of a sequence.

        Returns:
            Interval: A new Interval with reverse-complemented coordinates and strand.
        """
        if length is None:
            length = self.end
        return Interval(length - self.end, length - self.start, self.strand * -1)

    @classmethod
    def from_match(cls, item: Match, strand: Strand = Strand.UNSTRANDED) -> 'Interval':
        """Creates an Interval from a regular expression `Match` object.

        Args:
            item (Match): The regex match object, which has `start()` and `end()` methods.
            strand (Strand, optional): The strand to assign to the new interval. Defaults to UNSTRANDED.

        Returns:
            Interval: A new Interval corresponding to the match's span.
        """
        return cls(item.start(), item.end(), strand)

    @classmethod
    def from_int(cls, item: int, strand: Strand = Strand.UNSTRANDED, length: int | None = None) -> 'Interval':
        """Creates a 1-base-pair Interval from an integer coordinate.

        Args:
            item (int): The 0-based coordinate. If negative, it's treated as an offset from `length`.
            strand (Strand, optional): The strand to assign. Defaults to UNSTRANDED.
            length (int, optional): The parent sequence length, for resolving negative coordinates.

        Returns:
            Interval: A new Interval of length 1 at the specified position.
        """
        if item < 0 and length is not None: item += length
        return cls(item, item + 1, strand)

    @classmethod
    def from_slice(cls, item: slice, strand: Strand = Strand.UNSTRANDED, length: int | None = None) -> 'Interval':
        """Creates an Interval from a Python `slice` object.

        This allows for intuitive interval creation, e.g., `Interval.from_slice(10:20)`.

        Args:
            item (slice): The slice object with `start` and `stop` attributes.
            strand (Strand, optional): The strand to assign. Defaults to UNSTRANDED.
            length (int, optional): The parent sequence length, for resolving slices with `None` as the stop.

        Returns:
            Interval: A new Interval corresponding to the slice's range.

        Raises:
            ValueError: If the slice's `stop` is `None` and no `length` is provided.
        """
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
        """A universal factory method to create an Interval from various common types.

        This method acts as a type coercer, providing a single entry point for converting
        different representations (e.g., `slice`, `int`, `Match`) into a standardized `Interval` object.

        Args:
            item (IntervalLike): The object to convert. Can be an Interval, slice, int, or regex Match.
            strand (Strand, optional): The strand to assign if creating a new interval. Defaults to UNSTRANDED.
            length (int, optional): The parent sequence length, for resolving relative coordinates.

        Returns:
            Interval: The resulting Interval object.

        Raises:
            TypeError: If the input `item` is of an unsupported type.
        """
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
class Intervals:
    """A high-performance, vectorized SoA collection of genomic intervals.

    This class stores interval data in a Structure-of-Arrays (SoA) layout using NumPy arrays,
    which provides significant performance advantages over a list of individual `Interval` objects.
    By storing starts, ends, and strands in separate, contiguous arrays, this class enables:

    1.  **Memory Efficiency:** Reduced memory overhead compared to many small Python objects.
    2.  **Vectorized Operations:** NumPy's universal functions (ufuncs) can perform calculations
        (e.g., `expand`, `project`) on all intervals simultaneously without Python loops.
    3.  **Accelerated Queries:** Spatial queries are implemented with Numba-jitted kernels for C-like speed.

    The class is immutable by design. All methods that modify data (e.g., `sort`, `filter`, `expand`)
    return a new `Intervals` instance, leaving the original unchanged.

    Attributes:
        starts (npt.NDArray[np.int32]): A 1D array of 0-based start coordinates (inclusive).
        ends (npt.NDArray[np.int32]): A 1D array of 0-based end coordinates (exclusive).
        strands (npt.NDArray[np.int8]): A 1D array of strand identifiers (1, -1, or 0).
        original_indices (npt.NDArray[np.int32] | None): A 1D array tracking the original position of each
            interval before any sorting or filtering. This is crucial for relating query results back to
            their source data. If not provided, it's auto-generated as a simple range.
    """
    starts: npt.NDArray[np.int32]
    ends: npt.NDArray[np.int32]
    strands: npt.NDArray[np.int8]
    original_indices: npt.NDArray[np.int32] | None = None

    def __post_init__(self):
        """Initializes original_indices if they are not provided."""
        if self.original_indices is None:
            object.__setattr__(self, 'original_indices', np.arange(len(self.starts), dtype=np.int32))

    @classmethod
    def empty(cls) -> 'Intervals':
        """Creates an empty Intervals object with correctly typed, zero-length arrays.

        Returns:
            Intervals: An empty intervals instance.
        """
        return cls(
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8),
            np.empty(0, dtype=np.int32)
        )

    @classmethod
    def from_intervals(cls, intervals: Iterable[Interval]) -> 'Intervals':
        """Creates an Intervals object from an iterable of individual `Interval` objects.

        Args:
            intervals (Iterable[Interval]): A list, tuple, or generator of `Interval` objects.

        Returns:
            Intervals: A new collection containing the data from the intervals.
        """
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

    def __len__(self) -> int:
        """Returns the number of intervals in the batch."""
        return len(self.starts)

    def __dict__(self) -> dict:
        """Serializes the batch to a dictionary of lists."""
        return {'starts': self.starts.tolist(), 'ends': self.ends.tolist(), 'strands': self.strands.tolist()}

    @classmethod
    def from_dict(cls, d: dict) -> 'Intervals':
        """Deserializes intervals from a dictionary."""
        return cls(
            np.array(d['starts'], dtype=np.int32),
            np.array(d['ends'], dtype=np.int32),
            np.array(d['strands'], dtype=np.int8)
        )

    def __getitem__(self, item):
        """Accesses intervals by index, slice, or boolean mask.

        - If `item` is an integer, it returns a single `Interval` object for that position.
        - If `item` is a slice, list of indices, or a boolean mask, it returns a new,
          smaller `Intervals` collection containing only the selected intervals.

        Args:
            item: An integer index, slice, or boolean NumPy array.

        Returns:
            Interval | Intervals: A single interval or a new collection of intervals.

        Raises:
            IndexError: If an integer index is out of range.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return Interval(self.starts[item], self.ends[item], self.strands[item])

        return Intervals(
            self.starts[item],
            self.ends[item],
            self.strands[item],
            self.original_indices[item]
        )

    @classmethod
    def concat(cls, batches: Iterable['Intervals']) -> 'Intervals':
        """Concatenates multiple Intervals objects into a single, larger collection.

        Args:
            batches (Iterable[Intervals]): A list or iterable of collections to concatenate.

        Returns:
            Intervals: A new, combined intervals object.

        Raises:
            ValueError: If the input list of batches is empty.
        """
        batches = list(batches)
        if not batches:
            raise ValueError("Cannot concatenate empty list of batches")
        return cls(
            np.concatenate([b.starts for b in batches]),
            np.concatenate([b.ends for b in batches]),
            np.concatenate([b.strands for b in batches]),
            np.concatenate([b.original_indices for b in batches])
        )

    def shift(self, x: int | npt.NDArray[np.int32], y: int | npt.NDArray[np.int32] | None = None) -> 'Intervals':
        """Shifts all intervals in the collection by a given offset.
        
        Args:
            x (int | NDArray): The distance to shift the starts.
            y (int | NDArray, optional): The distance to shift the ends. If None, uses x.
            
        Returns:
            Intervals: A new Intervals collection with shifted coordinates.
        """
        if len(self) == 0:
            return self
            
        new_starts = self.starts + x
        new_ends = self.ends + (y if y is not None else x)
        
        return Intervals(
            np.asarray(new_starts, dtype=np.int32),
            np.asarray(new_ends, dtype=np.int32),
            self.strands,
            self.original_indices
        )

    def cluster_spatial(self, tolerance: int = 0,
                        group_by: npt.NDArray[np.integer] | None = None) -> npt.NDArray[np.int32]:
        """Clusters intervals that are within a spatial tolerance of one another.

        This performs a highly optimized 1D single-linkage clustering. It is ideal for
        grouping spatially close genes (e.g., into operons or BGCs) on a genome.

        Args:
            tolerance (int, optional): The maximum distance in base pairs between two intervals
                for them to be considered part of the same cluster. Defaults to 0.
            group_by (npt.NDArray[np.integer], optional): An integer array of the same length as
                the batch. Intervals will only be clustered together if they share the same
                group ID (e.g., the same contig). Defaults to None.

        Returns:
            npt.NDArray[np.int32]: A 1D array of cluster IDs parallel to the original intervals.
        """
        n = len(self)
        if n == 0:
            return np.empty(0, dtype=np.int32)

        if group_by is None:
            group_by = np.zeros(n, dtype=np.int32)
        else:
            group_by = np.asarray(group_by, dtype=np.int32)

        order = np.lexsort((self.ends, self.starts, group_by)).astype(np.int32)
        return _cluster_kernel(self.starts, self.ends, group_by, tolerance, order)

    def cluster_sequential(self, tolerance: int = 0, group_by: npt.NDArray[np.integer] | None = None,
                           enforce_strand: bool = False) -> npt.NDArray[np.int32]:
        """Clusters intervals based on their sequential index rather than physical distance.

        This performs 1D single-linkage clustering using the `original_indices` of the batch.
        It is ideal for grouping items like genes/ORFs by the number of intervening items between
        them, robustly bypassing large physical gaps (like IS elements).

        Args:
            tolerance (int, optional): The maximum number of allowed intervening intervals.
                A tolerance of 0 means intervals must be exactly adjacent. Defaults to 0.
            group_by (npt.NDArray[np.integer], optional): An integer array of the same length as
                the batch. Intervals will only be clustered together if they share the same
                group ID (e.g., the same contig). Defaults to None.
            enforce_strand (bool, optional): If True, intervals must also be on the same
                strand to be clustered together. Defaults to False.

        Returns:
            npt.NDArray[np.int32]: A 1D array of cluster IDs parallel to the original intervals.
        """
        n = len(self)
        if n == 0:
            return np.empty(0, dtype=np.int32)

        if group_by is None:
            group_by = np.zeros(n, dtype=np.int32)
        else:
            group_by = np.asarray(group_by, dtype=np.int32)

        if enforce_strand:
            order = np.lexsort((self.original_indices, self.strands, group_by)).astype(np.int32)
        else:
            order = np.lexsort((self.original_indices, group_by)).astype(np.int32)

        return _cluster_by_index_kernel(self.original_indices, group_by, self.strands, 
                                        tolerance, enforce_strand, order)

    def arrange(self, indices: npt.NDArray[np.integer],
                order: npt.NDArray[np.integer],
                starts: npt.NDArray[np.int32],
                ends: npt.NDArray[np.int32],
                strands: npt.NDArray[np.int8],
                gap: int = 500) -> 'Intervals':
        """Arranges interval coordinates across multiple disjoint pieces (e.g., contigs).
        
        Args:
            indices: Array mapping each interval to its piece index (0 to N-1).
            order: Array defining the layout order of the pieces.
            starts: Start coordinate of each piece.
            ends: End coordinate of each piece.
            strands: Orientation of each piece (1 for forward, -1 for reverse-complement).
            gap: Base pairs to insert between pieces.
            
        Returns:
            A new Intervals object with 1D arranged coordinates.
        """
        if len(self) == 0:
            return self

        n_pieces = len(starts)
        piece_plot_starts = np.zeros(n_pieces, dtype=np.int32)
        
        current_x = 0
        for i in order:
            p_len = ends[i] - starts[i]
            piece_plot_starts[i] = current_x
            current_x += p_len + gap

        new_starts = np.zeros(len(self), dtype=np.int32)
        new_ends = np.zeros(len(self), dtype=np.int32)
        new_strands = np.zeros(len(self), dtype=np.int8)

        # Vectorized coordinate transformation
        for i in range(n_pieces):
            mask = indices == i
            if not np.any(mask):
                continue
                
            p_s = starts[i]
            p_e = ends[i]
            orient = strands[i]
            offset = piece_plot_starts[i]
            
            g_s = self.starts[mask]
            g_e = self.ends[mask]
            g_str = self.strands[mask]
            
            if orient >= 0:
                new_starts[mask] = offset + (g_s - p_s)
                new_ends[mask] = offset + (g_e - p_s)
                new_strands[mask] = g_str
            else:
                new_starts[mask] = offset + (p_e - g_e)
                new_ends[mask] = offset + (p_e - g_s)
                new_strands[mask] = -g_str

        return Intervals(new_starts, new_ends, new_strands, self.original_indices)


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(cache=True, nogil=True)
def _cluster_kernel(starts: npt.NDArray[np.int32], ends: npt.NDArray[np.int32],
                    groups: npt.NDArray[np.int32], tolerance: int,
                    order: npt.NDArray[np.int32]) -> npt.NDArray[np.int32]:
    """Numba-accelerated 1D spatial clustering."""
    n = len(starts)
    cluster_ids = np.empty(n, dtype=np.int32)
    if n == 0: return cluster_ids

    curr_cluster = 0
    first_idx = order[0]
    curr_e = ends[first_idx]
    curr_g = groups[first_idx]
    cluster_ids[first_idx] = curr_cluster

    for i in range(1, n):
        idx = order[i]
        s, e, g = starts[idx], ends[idx], groups[idx]

        if g == curr_g and s <= curr_e + tolerance:
            curr_e = max(curr_e, e)
        else:
            curr_cluster += 1
            curr_e = e
            curr_g = g

        cluster_ids[idx] = curr_cluster

    return cluster_ids


@njit(cache=True, nogil=True)
def _cluster_by_index_kernel(indices: npt.NDArray[np.int32], groups: npt.NDArray[np.int32],
                             strands: npt.NDArray[np.int8], tolerance: int,
                             enforce_strand: bool, order: npt.NDArray[np.int32]) -> npt.NDArray[np.int32]:
    """Numba-accelerated 1D index-based clustering."""
    n = len(indices)
    cluster_ids = np.empty(n, dtype=np.int32)
    if n == 0: return cluster_ids

    curr_cluster = 0
    first_idx = order[0]
    curr_i = indices[first_idx]
    curr_g = groups[first_idx]
    curr_st = strands[first_idx]
    cluster_ids[first_idx] = curr_cluster

    for i in range(1, n):
        idx = order[i]
        val_i = indices[idx]
        g = groups[idx]
        st = strands[idx]

        same_group = (g == curr_g)
        same_strand = (st == curr_st) if enforce_strand else True

        if same_group and same_strand and (abs(val_i - curr_i) <= tolerance + 1):
            curr_i = max(curr_i, val_i)
        else:
            curr_cluster += 1
            curr_i = val_i
            curr_g = g
            curr_st = st

        cluster_ids[idx] = curr_cluster

    return cluster_ids
