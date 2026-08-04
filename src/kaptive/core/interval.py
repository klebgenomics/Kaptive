r"""Genomic interval representation, strand orientations, and vectorized Structure-of-Arrays operations.

This module provides discrete coordinate intervals via [`Interval`][kaptive.core.interval.Interval],
strand representations via [`Strand`][kaptive.core.interval.Strand], and high-performance,
vectorized interval collections via [`Intervals`][kaptive.core.interval.Intervals].
Includes JIT-compiled Numba kernels for 1D spatial/sequential clustering and overlap culling.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import IntEnum
from re import Match
from typing import Any, Self

import numpy as np
import numpy.typing as npt
from numba import njit

from kaptive.core.collections import BatchedContainer


# Enums ----------------------------------------------------------------------------------------------------------------
class Strand(IntEnum):
    r"""Integer representation of genomic strand orientation.

    Supports conversion from common string formats (`'+'`, `'-'`, `'1'`, `'-1'`).

    Attributes:
        FORWARD (int): Forward strand (`+1`).
        REVERSE (int): Reverse strand (`-1`).
        UNSTRANDED (int): Unstranded or unknown orientation (`0`).
    """

    FORWARD = 1
    REVERSE = -1
    UNSTRANDED = 0

    @classmethod
    def _missing_(cls, value: object) -> Strand:
        r"""Coerce strings, bytes, or integers to a valid [`Strand`][kaptive.core.interval.Strand] instance.

        Args:
            value (object): Value to coerce (e.g. `'+'`, `'-'`, `b'+'`, `1`, `-1`).

        Returns:
            Strand: Corresponding [`Strand`][kaptive.core.interval.Strand] enum member.
        """
        if isinstance(value, bytes):
            value = value.decode("ascii")
        if isinstance(value, str):
            if value in ("+", "1", "+1"):
                return Strand.FORWARD
            elif value in ("-", "-1"):
                return Strand.REVERSE
        return Strand.UNSTRANDED

    def __str__(self) -> str:
        r"""Return the string representation of the strand (`'+'`, `'-'`, or `'.'`).

        Returns:
            str: Single-character string representation.
        """
        if self == Strand.FORWARD:
            return "+"
        if self == Strand.REVERSE:
            return "-"
        return "."


# Classes --------------------------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Interval:
    r"""A single 0-based, half-open genomic interval with strand orientation.

    Represents a discrete sequence segment defined by `start` (inclusive), `end` (exclusive),
    and `strand` orientation.

    Attributes:
        start (int): 0-based starting position (inclusive).
        end (int): 0-based ending position (exclusive).
        strand (Strand): Strand orientation [`Strand`][kaptive.core.interval.Strand].
    """

    start: int
    end: int
    strand: Strand = Strand.UNSTRANDED

    def __len__(self) -> int:
        r"""Calculate the length of the interval in base pairs (`end - start`).

        Returns:
            int: Interval length.
        """
        return self.end - self.start

    def __contains__(self, item: IntervalLike) -> bool:
        r"""Check if a coordinate or another interval is fully contained within this interval.

        Args:
            item (IntervalLike): Integer coordinate, slice, regex Match, or
                [`Interval`][kaptive.core.interval.Interval].

        Returns:
            bool: True if item is completely bounded within start and end coordinates.
        """
        if isinstance(item, int):
            return self.start <= item < self.end
        interval_obj = Interval.from_item(item)
        return self.start <= interval_obj.start and self.end >= interval_obj.end

    def __add__(self, other: IntervalLike) -> Interval:
        r"""Compute the minimal bounding interval covering both this interval and another.

        Args:
            other (IntervalLike): Interval-like object to merge with.

        Returns:
            Interval: A new [`Interval`][kaptive.core.interval.Interval] representing the combined span.
        """
        other_obj = Interval.from_item(other)
        new_strand = self.strand if self.strand == other_obj.strand else Strand.UNSTRANDED
        return Interval(min(self.start, other_obj.start), max(self.end, other_obj.end), new_strand)

    def __radd__(self, other: IntervalLike) -> Interval:
        r"""Reverse addition operator supporting `other + self`.

        Args:
            other (IntervalLike): Interval-like object to merge with.

        Returns:
            Interval: Combined bounding [`Interval`][kaptive.core.interval.Interval].
        """
        return self.__add__(other)

    def shift(self, x: int, y: int | None = None) -> Interval:
        r"""Shift interval coordinates by fixed distances.

        Args:
            x (int): Distance to shift the start coordinate.
            y (int | None): Distance to shift the end coordinate. If None, defaults to `x`.

        Returns:
            Interval: A new shifted [`Interval`][kaptive.core.interval.Interval].
        """
        return Interval(self.start + x, self.end + (y if y is not None else x), self.strand)

    def expand(self, left: int, right: int, clip_length: int | None = None) -> Interval:
        r"""Expand or shrink the interval by specified base-pair amounts.

        Args:
            left (int): Amount to extend start leftward (subtract from `start`).
            right (int): Amount to extend end rightward (add to `end`).
            clip_length (int | None): Maximum boundary constraint length.

        Returns:
            Interval: A new expanded [`Interval`][kaptive.core.interval.Interval].
        """
        new_start = max(0, self.start - left)
        new_end = self.end + right
        if clip_length is not None:
            new_end = min(new_end, clip_length)

        return Interval(new_start, new_end, self.strand)

    def reverse_complement(self, length: int | None = None) -> Interval:
        r"""Calculate the mirrored coordinates of the interval on the opposite strand.

        Args:
            length (int | None): Total length of parent sequence. If None, defaults to `end`.

        Returns:
            Interval: Mirrored [`Interval`][kaptive.core.interval.Interval] on opposite strand.
        """
        if length is None:
            length = self.end
        return Interval(length - self.end, length - self.start, Strand(self.strand * -1))

    @classmethod
    def from_match(cls, item: Match, strand: Strand = Strand.UNSTRANDED) -> Interval:  # type: ignore
        r"""Create an [`Interval`][kaptive.core.interval.Interval] from a regular expression match object.

        Args:
            item (Match): Regex match object containing `start()` and `end()`.
            strand (Strand): Strand orientation [`Strand`][kaptive.core.interval.Strand].

        Returns:
            Interval: A new [`Interval`][kaptive.core.interval.Interval].
        """
        return cls(item.start(), item.end(), strand)

    @classmethod
    def from_int(cls, item: int, strand: Strand = Strand.UNSTRANDED, length: int | None = None) -> Interval:
        r"""Create a 1-bp [`Interval`][kaptive.core.interval.Interval] from an integer coordinate.

        Args:
            item (int): 0-based coordinate.
            strand (Strand): Strand orientation [`Strand`][kaptive.core.interval.Strand].
            length (int | None): Parent sequence length to resolve negative indices.

        Returns:
            Interval: 1-bp [`Interval`][kaptive.core.interval.Interval].
        """
        if item < 0 and length is not None:
            item += length
        return cls(item, item + 1, strand)

    @classmethod
    def from_slice(cls, item: slice, strand: Strand = Strand.UNSTRANDED, length: int | None = None) -> Interval:
        r"""Create an [`Interval`][kaptive.core.interval.Interval] from a Python slice object.

        Args:
            item (slice): Python slice with `start` and `stop` values.
            strand (Strand): Strand orientation [`Strand`][kaptive.core.interval.Strand].
            length (int | None): Parent sequence length for open-ended slices.

        Returns:
            Interval: A new [`Interval`][kaptive.core.interval.Interval].

        Raises:
            ValueError: If slice `stop` is None and no `length` parameter is provided.
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
    def from_item(cls, item: IntervalLike, strand: Strand = Strand.UNSTRANDED, length: int | None = None) -> Interval:
        r"""Universal factory coercing various representations into an [`Interval`][kaptive.core.interval.Interval].

        Args:
            item (IntervalLike): Object to convert (Interval, slice, int, or regex Match).
            strand (Strand): Strand orientation [`Strand`][kaptive.core.interval.Strand].
            length (int | None): Parent sequence length for coordinate resolution.

        Returns:
            Interval: Standardized [`Interval`][kaptive.core.interval.Interval] instance.

        Raises:
            TypeError: If `item` is of an unsupported type.
        """
        if isinstance(item, cls):
            return item
        if (interval := getattr(item, "interval", None)) is not None:
            return interval
        if isinstance(item, Match):
            return cls.from_match(item, strand)
        if isinstance(item, int):
            return cls.from_int(item, strand, length)
        if isinstance(item, slice):
            return cls.from_slice(item, strand, length)
        raise TypeError(item)


IntervalLike = slice | int | Match | Interval


@dataclass(frozen=True, slots=True)
class Intervals(BatchedContainer[Interval, "Intervals"]):
    r"""A high-performance Structure-of-Arrays (SoA) collection of genomic intervals.

    Stores starts, ends, strands, and original tracking indices as parallel 1D NumPy arrays
    for memory-efficient vectorized operations and Numba spatial algorithms.

    Attributes:
        starts (npt.NDArray[np.int32]): 1D array of 0-based start coordinates.
        ends (npt.NDArray[np.int32]): 1D array of 0-based end coordinates.
        strands (npt.NDArray[np.int8]): 1D array of strand values (+1, -1, 0).
        original_indices (npt.NDArray[np.int32] | None): 1D array tracking original item indices.
    """

    starts: npt.NDArray[np.int32]
    ends: npt.NDArray[np.int32]
    strands: npt.NDArray[np.int8]
    original_indices: npt.NDArray[np.int32] | None = None

    def __post_init__(self) -> None:
        r"""Initialize default `original_indices` array if not explicitly supplied."""
        if self.original_indices is None:
            object.__setattr__(self, "original_indices", np.arange(len(self.starts), dtype=np.int32))

    @classmethod
    def empty(cls) -> Intervals:
        r"""Create an empty [`Intervals`][kaptive.core.interval.Intervals] object with zero-length arrays.

        Returns:
            Intervals: Empty [`Intervals`][kaptive.core.interval.Intervals] collection.
        """
        return cls(
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8),
            np.empty(0, dtype=np.int32),
        )

    @classmethod
    def from_intervals(cls, intervals: Iterable[Interval]) -> Intervals:
        r"""Construct an [`Intervals`][kaptive.core.interval.Intervals] collection from individual intervals.

        Args:
            intervals (Iterable[Interval]): Iterable of [`Interval`][kaptive.core.interval.Interval] instances.

        Returns:
            Intervals: Vectorized [`Intervals`][kaptive.core.interval.Intervals] collection.
        """
        # OPTIMIZATION: Fast C-level list comprehension + zip extraction
        data = [(i.start, i.end, i.strand) for i in intervals]
        if not data:
            return cls.empty()

        start_vals, end_vals, strand_vals = zip(*data, strict=False)
        return cls(
            np.array(start_vals, dtype=np.int32),
            np.array(end_vals, dtype=np.int32),
            np.array(strand_vals, dtype=np.int8),
        )

    def __len__(self) -> int:
        r"""Return the number of intervals in the batch.

        Returns:
            int: Total interval count.
        """
        return len(self.starts)

    def to_dict(self) -> dict[str, list]:  # type: ignore
        r"""Serialize the interval batch into a dictionary of python lists.

        Returns:
            dict[str, list]: Dictionary containing `'starts'`, `'ends'`, and `'strands'`.
        """
        return {"starts": self.starts.tolist(), "ends": self.ends.tolist(), "strands": self.strands.tolist()}

    @classmethod
    def from_dict(cls, d: dict) -> Intervals:  # type: ignore
        r"""Deserialize an [`Intervals`][kaptive.core.interval.Intervals] collection from a dictionary.

        Args:
            d (dict): Dictionary with `'starts'`, `'ends'`, and `'strands'` keys.

        Returns:
            Intervals: Restored [`Intervals`][kaptive.core.interval.Intervals] collection.
        """
        return cls(
            np.array(d["starts"], dtype=np.int32),
            np.array(d["ends"], dtype=np.int32),
            np.array(d["strands"], dtype=np.int8),
        )

    def __getitem__(self, item: Any) -> Interval | Intervals:
        r"""Access intervals by integer index, slice, or boolean/integer NumPy array.

        Args:
            item (Any): Integer index, slice, or NumPy mask array.

        Returns:
            Interval | Intervals: A single [`Interval`][kaptive.core.interval.Interval] or a
                sliced [`Intervals`][kaptive.core.interval.Intervals] collection.

        Raises:
            IndexError: If integer index is out of bounds.
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
            self.original_indices[item] if self.original_indices is not None else None,
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenate multiple [`Intervals`][kaptive.core.interval.Intervals] collections into a single batch.

        Args:
            batches (Iterable[Intervals]): Iterable of [`Intervals`][kaptive.core.interval.Intervals] objects to merge.

        Returns:
            Intervals: Concatenated [`Intervals`][kaptive.core.interval.Intervals] collection.

        Raises:
            ValueError: If `batches` is empty.
        """
        batches_list = list(batches)
        if not batches_list:
            raise ValueError("Cannot concatenate empty list of batches")
        return cls(
            np.concatenate([b.starts for b in batches_list]),
            np.concatenate([b.ends for b in batches_list]),
            np.concatenate([b.strands for b in batches_list]),
            np.concatenate([b.original_indices for b in batches_list])
            if batches_list[0].original_indices is not None
            else None,
        )

    def shift(self, x: int | npt.NDArray[np.int32], y: int | npt.NDArray[np.int32] | None = None) -> Intervals:
        r"""Vectorized coordinate shift of all intervals in the collection.

        Args:
            x (int | npt.NDArray[np.int32]): Shift distance for starts.
            y (int | npt.NDArray[np.int32] | None): Shift distance for ends. Defaults to `x`.

        Returns:
            Intervals: Shifted [`Intervals`][kaptive.core.interval.Intervals] collection.
        """
        if len(self) == 0:
            return self

        new_starts = self.starts + x
        new_ends = self.ends + (y if y is not None else x)

        return Intervals(
            np.asarray(new_starts, dtype=np.int32),
            np.asarray(new_ends, dtype=np.int32),
            self.strands,
            self.original_indices,
        )

    def cull_overlaps(
        self,
        order: npt.NDArray[np.int32],
        max_overlap_fraction: float = 0.1,
        group_by: npt.NDArray[np.integer] | None = None,
        secondary_group_by: npt.NDArray[np.integer] | None = None,
    ) -> npt.NDArray[np.bool_]:
        r"""Greedily cull overlapping intervals based on prior ordering and max overlap thresholds.

        Args:
            order (npt.NDArray[np.int32]): Priority evaluation order indices.
            max_overlap_fraction (float): Maximum allowable overlap ratio. Defaults to 0.1.
            group_by (npt.NDArray[np.integer] | None): Primary grouping array (e.g. contig ID).
            secondary_group_by (npt.NDArray[np.integer] | None): Secondary grouping array.

        Returns:
            npt.NDArray[np.bool_]: Boolean mask array indicating kept intervals.
        """
        n = len(self)
        if n == 0:
            return np.empty(0, dtype=np.bool_)

        if group_by is None:
            group_by_arr = np.zeros(n, dtype=np.int32)
        else:
            group_by_arr = np.asarray(group_by, dtype=np.int32)

        if secondary_group_by is None:
            secondary_group_by_arr = np.zeros(n, dtype=np.int32)
        else:
            secondary_group_by_arr = np.asarray(secondary_group_by, dtype=np.int32)

        return _cull_overlaps_kernel(
            order, group_by_arr, secondary_group_by_arr, self.starts, self.ends, max_overlap_fraction, n
        )

    def cluster_spatial(
        self, tolerance: int = 0, group_by: npt.NDArray[np.integer] | None = None
    ) -> npt.NDArray[np.int32]:
        r"""Perform fast 1D single-linkage spatial clustering of intervals.

        Args:
            tolerance (int): Maximum base-pair distance threshold for clustering. Defaults to 0.
            group_by (npt.NDArray[np.integer] | None): Grouping array (e.g. contig ID).

        Returns:
            npt.NDArray[np.int32]: 1D array of cluster IDs parallel to original intervals.
        """
        n = len(self)
        if n == 0:
            return np.empty(0, dtype=np.int32)

        if group_by is None:
            group_by_arr = np.zeros(n, dtype=np.int32)
        else:
            group_by_arr = np.asarray(group_by, dtype=np.int32)

        order = np.lexsort((self.ends, self.starts, group_by_arr)).astype(np.int32)
        return _cluster_kernel(self.starts, self.ends, group_by_arr, tolerance, order)

    def cluster_sequential(
        self,
        tolerance: int = 0,
        group_by: npt.NDArray[np.integer] | None = None,
        enforce_strand: bool = False,
    ) -> npt.NDArray[np.int32]:
        r"""Perform index-based sequential clustering independent of physical base-pair distance.

        Args:
            tolerance (int): Maximum allowed intervening item count. Defaults to 0.
            group_by (npt.NDArray[np.integer] | None): Grouping array (e.g. contig ID).
            enforce_strand (bool): If True, restricts clustering to identical strands.

        Returns:
            npt.NDArray[np.int32]: 1D array of cluster IDs parallel to original intervals.
        """
        n = len(self)
        if n == 0:
            return np.empty(0, dtype=np.int32)

        if group_by is None:
            group_by_arr = np.zeros(n, dtype=np.int32)
        else:
            group_by_arr = np.asarray(group_by, dtype=np.int32)

        indices = self.original_indices if self.original_indices is not None else np.zeros(n, dtype=np.int32)

        if enforce_strand:
            order = np.lexsort((indices, self.strands, group_by_arr)).astype(np.int32)
        else:
            order = np.lexsort((indices, group_by_arr)).astype(np.int32)

        return _cluster_by_index_kernel(indices, group_by_arr, self.strands, tolerance, enforce_strand, order)

    def arrange(
        self,
        indices: npt.NDArray[np.integer],
        order: npt.NDArray[np.integer],
        starts: npt.NDArray[np.int32],
        ends: npt.NDArray[np.int32],
        strands: npt.NDArray[np.int8],
        gap: int = 500,
    ) -> Intervals:
        r"""Arrange interval coordinates across disjoint contig/piece layouts into a 1D space.

        Args:
            indices (npt.NDArray[np.integer]): Piece index mapping for each interval.
            order (npt.NDArray[np.integer]): Piece layout order.
            starts (npt.NDArray[np.int32]): Piece start coordinates.
            ends (npt.NDArray[np.int32]): Piece end coordinates.
            strands (npt.NDArray[np.int8]): Piece orientations (+1 or -1).
            gap (int): Base-pair padding between adjacent pieces. Defaults to 500.

        Returns:
            Intervals: Transformed 1D arranged [`Intervals`][kaptive.core.interval.Intervals] collection.
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
            orient = Strand(strands[i])
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
def _cluster_kernel(
    starts: npt.NDArray[np.int32],
    ends: npt.NDArray[np.int32],
    groups: npt.NDArray[np.int32],
    tolerance: int,
    order: npt.NDArray[np.int32],
) -> npt.NDArray[np.int32]:
    r"""Numba-accelerated 1D spatial clustering.

    Args:
        starts (npt.NDArray[np.int32]): Start coordinates array.
        ends (npt.NDArray[np.int32]): End coordinates array.
        groups (npt.NDArray[np.int32]): Grouping identifiers array.
        tolerance (int): Maximum base pair distance allowed.
        order (npt.NDArray[np.int32]): Sorted order index array.

    Returns:
        npt.NDArray[np.int32]: Cluster IDs parallel to original intervals.
    """
    n = len(starts)
    cluster_ids = np.empty(n, dtype=np.int32)
    if n == 0:
        return cluster_ids

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
def _cluster_by_index_kernel(
    indices: npt.NDArray[np.int32],
    groups: npt.NDArray[np.int32],
    strands: npt.NDArray[np.int8],
    tolerance: int,
    enforce_strand: bool,
    order: npt.NDArray[np.int32],
) -> npt.NDArray[np.int32]:
    r"""Numba-accelerated 1D index-based clustering.

    Args:
        indices (npt.NDArray[np.int32]): Original item tracking indices.
        groups (npt.NDArray[np.int32]): Grouping identifiers array.
        strands (npt.NDArray[np.int8]): Strand orientations (+1, -1, 0).
        tolerance (int): Maximum allowed intervening item count.
        enforce_strand (bool): Whether to enforce identical strands.
        order (npt.NDArray[np.int32]): Sorted order index array.

    Returns:
        npt.NDArray[np.int32]: Cluster IDs parallel to original intervals.
    """
    n = len(indices)
    cluster_ids = np.empty(n, dtype=np.int32)
    if n == 0:
        return cluster_ids

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

        same_group = g == curr_g
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


@njit(cache=True, nogil=True)
def _cull_overlaps_kernel(
    order: npt.NDArray[np.int32],
    group1: npt.NDArray[np.int32],
    group2: npt.NDArray[np.int32],
    starts: npt.NDArray[np.int32],
    ends: npt.NDArray[np.int32],
    max_overlap_fraction: float,
    n: int,
) -> npt.NDArray[np.bool_]:
    r"""Numba-accelerated greedy overlap culling.

    Args:
        order (npt.NDArray[np.int32]): Evaluation order indices.
        group1 (npt.NDArray[np.int32]): Primary group IDs.
        group2 (npt.NDArray[np.int32]): Secondary group IDs.
        starts (npt.NDArray[np.int32]): Start coordinates.
        ends (npt.NDArray[np.int32]): End coordinates.
        max_overlap_fraction (float): Maximum allowed overlap ratio.
        n (int): Number of intervals.

    Returns:
        npt.NDArray[np.bool_]: Boolean mask array of kept intervals.
    """
    kept_mask = np.zeros(n, dtype=np.bool_)

    for i in range(n):
        idx = order[i]
        g1 = group1[idx]
        g2 = group2[idx]
        s = starts[idx]
        e = ends[idx]
        length = e - s

        if length <= 0:
            continue

        overlap_found = False
        # Check against previously kept intervals in O(N^2) (Very fast in pure C)
        for j in range(i):
            prev_idx = order[j]
            if not kept_mask[prev_idx] or group1[prev_idx] != g1 or group2[prev_idx] != g2:
                continue

            ks, ke = starts[prev_idx], ends[prev_idx]
            overlap = min(e, ke) - max(s, ks)
            if overlap > 0 and (overlap / min(length, ke - ks)) > max_overlap_fraction:
                overlap_found = True
                break

        if not overlap_found:
            kept_mask[idx] = True

    return kept_mask
