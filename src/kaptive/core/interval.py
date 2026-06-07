"""
Genomic interval representation with strand and context, plus batched interval operations.
"""
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Self, Union, Match
from enum import IntEnum, auto

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
        if self == Strand.FORWARD: return '+'
        if self == Strand.REVERSE: return '-'
        return '.'


class Context(IntEnum):
    """
    Spatial relationship between two genomic intervals.

    Used to describe the position of a 'passenger' or 'flank' gene relative
    to a primary target (e.g. a mobile element insertion site).
    """
    UPSTREAM = auto()
    DOWNSTREAM = auto()
    INSIDE = auto()
    OVERLAPPING = auto()
    OVERLAPPING_START = auto()
    OVERLAPPING_END = auto()


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

    def relate(self, other: IntervalLike) -> Context:
        """Calculates the spatial relationship between this interval and another.

        This method determines how two intervals are positioned relative to each other, considering strand.
        For example, it can identify if another interval is upstream, downstream, overlapping, or contained within
        this one.

        Args:
            other (IntervalLike): The other interval to compare against.

        Returns:
            Context: A `Context` enum value (e.g., UPSTREAM, INSIDE, OVERLAPPING_5_PRIME) describing
                the relationship.
        """
        other = Interval.from_item(other)
        return Context(_core_relate(self.start, self.end, self.strand, other.start, other.end))

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
    def empty(cls) -> Self:
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
    def from_intervals(cls, intervals: Iterable[Interval]) -> Self:
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

    def max_len(self) -> int:
        """Returns the length of the longest interval in the batch.

        Returns:
            int: The maximum length, or 0 if the batch is empty.
        """
        return np.max(self.ends - self.starts) if len(self) > 0 else 0

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

    def copy(self):
        """Creates a deep copy, including all underlying NumPy arrays.

        Returns:
            Intervals: A new `Intervals` object with copied data.
        """
        return Intervals(
            self.starts.copy(),
            self.ends.copy(),
            self.strands.copy(),
            self.original_indices.copy()
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:
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

    def sort(self) -> 'Intervals':
        """Returns a new collection sorted by start coordinate, then end coordinate.

        This is a prerequisite for many efficient spatial algorithms (e.g., merging, sweeping).
        The implementation includes a fast Numba-jitted check to skip sorting if the batch is
        already sorted.

        Returns:
            Intervals: A new, sorted `Intervals` collection.
        """
        if len(self) < 2 or _is_sorted_kernel(self.starts, self.ends):
            return self

        order = np.lexsort((self.ends, self.starts))
        return Intervals(
            self.starts[order],
            self.ends[order],
            self.strands[order],
            self.original_indices[order]
        )

    def filter(self, mask: IntervalLike) -> 'Intervals':
        """Returns a new collection containing only intervals matched by the mask.

        This is a convenience wrapper around `__getitem__` for clarity.

        Args:
            mask (IntervalLike): A slice, list of indices, or boolean array.

        Returns:
            Intervals: A new, filtered collection.
        """
        if isinstance(mask, (slice, int, np.integer)):
            if isinstance(mask, (int, np.integer)):
                mask = [mask]
            return self[mask]
        return self[np.asarray(mask)]

    @property
    def envelope(self) -> Interval | None:
        """Returns a single bounding `Interval` encompassing the absolute min and max coordinates of the batch.
        
        Returns None if the batch is empty.
        """
        if len(self) == 0:
            return None
        
        # Assumes all intervals are on the same sequence/contig
        return Interval(int(np.min(self.starts)), int(np.max(self.ends)))

    @property
    def centers(self) -> npt.NDArray:
        """Vectorized calculation of all interval midpoints.

        Returns:
            npt.NDArray: A float array of the midpoints.
        """
        return (self.starts + self.ends) / 2

    @property
    def lengths(self) -> npt.NDArray[np.int32]:
        """Vectorized calculation of all interval lengths.

        Returns:
            npt.NDArray[np.int32]: An integer array of the lengths.
        """
        return self.ends - self.starts

    def query(self, start: int, end: int) -> npt.NDArray[np.int32]:
        """Finds indices of all intervals in the batch that overlap with a given query range.

        This method uses a highly optimized, Numba-jitted kernel. It first performs a binary
        search (`searchsorted`) to narrow down the potential candidates and then iterates backward
        to check for overlaps. This is significantly faster than a linear scan for large batches.
        The batch should ideally be sorted for best performance.

        Args:
            start (int): The start of the query range.
            end (int): The end of the query range.

        Returns:
            np.ndarray: A 1D integer array of the indices of overlapping intervals.
        """
        if len(self) == 0:
            return np.empty(0, dtype=np.int32)
        return _query_kernel(self.starts, self.ends, start, end, self.max_len())

    def merge(self, tolerance: int = 0) -> 'Intervals':
        """Merges overlapping or adjacent intervals into single bounding boxes.

        This method collapses all intervals that overlap or are within `tolerance` base pairs of
        each other into a single, minimal bounding interval. The batch must be sorted first.
        The operation is performed by a fast Numba-jitted kernel.

        Args:
            tolerance (int, optional): The maximum gap between intervals to be considered
                for merging. Defaults to 0.

        Returns:
            Intervals: A new collection with the merged intervals. Note that `original_indices`
                are not preserved in the merged output as they become ambiguous.
        """
        if len(self) == 0:
            return self
        out = _merge_kernel(self.starts, self.ends, self.strands, tolerance)
        return type(self)(out[0], out[1], out[2])

    def expand(self, left: int | npt.NDArray, right: int | npt.NDArray, 
               clip_lengths: int | npt.NDArray | None = None) -> 'Intervals':
        """Vectorized expansion (or shrinkage) of all intervals in the collection.

        This method applies an expansion to every interval simultaneously. The `left` and `right`
        parameters can be single integers (applied to all) or NumPy arrays of the same length
        as the batch for per-interval adjustments.

        Args:
            left (int | npt.NDArray): Amount to expand the start of each interval.
            right (int | npt.NDArray): Amount to expand the end of each interval.
            clip_lengths (int | npt.NDArray | None, optional): Maximum length(s) to clip the
                end coordinates against. Can be a single int or an array. Defaults to None.

        Returns:
            Intervals: A new collection with the expanded intervals.
        """
        new_starts = np.maximum(0, self.starts - left)
        new_ends = self.ends + right
        
        if clip_lengths is not None:
            new_ends = np.minimum(new_ends, clip_lengths)
            
        return Intervals(
            starts=new_starts.astype(self.starts.dtype),
            ends=new_ends.astype(self.ends.dtype),
            strands=self.strands.copy(),
            original_indices=self.original_indices.copy()
        )

    def project(self, shift: int, flip_length: int | None = None) -> 'Intervals':
        """Applies a coordinate transformation to all intervals in the collection.

        This is a vectorized version of `Interval.shift` and `Interval.reverse_complement`. It can
        shift all intervals and optionally mirror them within a given length, which is useful for
        projecting coordinates between different reference frames (e.g., from a gene to a contig,
        or from a contig's forward strand to its reverse strand).

        Args:
            shift (int): The distance to add to all start and end coordinates.
            flip_length (int, optional): If provided, coordinates are mirrored within this length
                (e.g., `new_start = flip_length - old_end`) and strands are inverted. This is
                used for reverse-complementation. Defaults to None.

        Returns:
            Intervals: A new collection with the transformed coordinates.
        """
        if flip_length is not None:
            new_starts = (flip_length - self.ends).astype(self.starts.dtype)
            new_ends = (flip_length - self.starts).astype(self.ends.dtype)
            new_strands = (-self.strands).astype(self.strands.dtype)
        else:
            new_starts = self.starts.copy()
            new_ends = self.ends.copy()
            new_strands = self.strands.copy()

        new_starts += shift
        new_ends += shift

        return Intervals(
            starts=new_starts,
            ends=new_ends,
            strands=new_strands,
            original_indices=self.original_indices.copy()
        )

    def overlaps_with(self, other: 'Intervals') -> np.ndarray:
        """Performs a vectorized, many-to-many overlap check between this collection and another.

        This method constructs a 2D boolean matrix where `matrix[i, j]` is true if `self[i]`
        overlaps with `other[j]`. It then reduces this matrix along the `other` axis (`axis=1`)
        to produce a 1D boolean mask indicating which intervals in `self` overlap with *any*
        interval in `other`. This is highly efficient for checking overlaps between two large sets.

        Args:
            other (Intervals): The other collection to check for overlaps against.

        Returns:
            np.ndarray: A 1D boolean array of the same length as `self`, where `True` at index `i`
                means `self[i]` overlaps with at least one interval in `other`.
        """
        if len(self) == 0 or len(other) == 0:
            return np.zeros(len(self), dtype=bool)

        overlaps = (self.starts[:, None] < other.ends[None, :]) & \
                   (other.starts[None, :] < self.ends[:, None])
        return overlaps.any(axis=1)

    def cluster(self, tolerance: int = 0, group_by: npt.NDArray[np.integer] | None = None) -> npt.NDArray[np.int32]:
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

    def cluster_by_index(self, tolerance: int = 0, group_by: npt.NDArray[np.integer] | None = None,
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
