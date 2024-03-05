"""
Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/klebgenomics/Kaptive

This file is part of Kaptive. Kaptive is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kaptive is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kaptive.
If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import annotations

from typing import Iterable, Generator


class IntRangeError(Exception):
    pass


class IntRange(object):
    """
    Modified from: https://github.com/rrwick/Verticall/blob/main/verticall/intrange.py
    This class contains one or more integer ranges in the form of sorted integer tuples.
    Overlapping ranges will be merged together and stored in a Python-like fashion where the last value in each range is
    exclusive.
    """

    def __init__(self, ranges: Iterable[tuple[int, int]] | tuple[int, int] | None = None):
        self._ranges = []
        if not ranges:
            pass
        elif isinstance(ranges, (tuple, Iterable)):
            self.add_ranges([ranges])
        else:
            raise TypeError(f'Cannot create IntRange from {type(ranges)}')

    def add_range(self, start: int, end: int):
        """Creates an IntRange object from a single range."""
        if not isinstance(start, int) or not isinstance(end, int):
            TypeError('Start and end must be integers')
        self._ranges = merge_ranges(self._ranges + [(start, end)] if start < end else [(end, start)])

    def add_ranges(self, ranges: Iterable[tuple[int, int]]):
        """Creates an IntRange object from a list of ranges."""
        if not isinstance(ranges, Iterable):
            raise TypeError('Ranges must be iterable')
        for r in ranges:
            if not isinstance(r, tuple) and len(r) != 2:
                raise TypeError('Ranges must be tuples of length 2')
            self.add_range(*r)

    def __repr__(self):
        return str(self._ranges)

    def __len__(self):
        return sum(x[1] - x[0] for x in self._ranges) if self._ranges else 0

    def __iter__(self):
        return iter(self._ranges)

    def __add__(self, other: IntRange | tuple[int, int] | list[tuple[int, int]]):
        if isinstance(other, IntRange):
            return IntRange(self._ranges + other._ranges)
        elif isinstance(other, tuple):
            return IntRange(self._ranges + [other])
        elif isinstance(other, list):
            return IntRange(self._ranges + other)

    def __contains__(self, item: int | float | tuple[int, int] | list[tuple[int, int]] | IntRange) -> bool:
        if not isinstance(item, (int, float, tuple, list, IntRange)):
            raise TypeError(f'Cannot check if {type(item)} is in IntRange')
        if isinstance(item, (int, float)):
            return any([i[0] <= item <= i[1] for i in self._ranges])
        item = IntRange(item) if not isinstance(item, IntRange) else item
        for r in item._ranges:
            for x in self._ranges:
                if x[0] <= r[0] <= x[1] and x[0] <= r[1] <= x[1]:
                    return True
        return False

    def overlaps(self, other: IntRange) -> bool:
        """
        Returns whether this IntRange overlaps with another IntRange using the range_overlap function.
        """
        for other_range in other._ranges:
            for self_range in self._ranges:
                if range_overlap(self_range, other_range) > 0:
                    return True
        return False


def merge_ranges(ranges: list[tuple[int | float, int | float]], tolerance: int | float = 0, skip_sort: bool = False) -> Generator[tuple[int | float, int | float]]:
    """
    Merge overlapping ranges
    :param ranges: List of tuples of start and end positions
    :param tolerance: Integer or float of tolerance for merging ranges
    :param skip_sort: Skip sorting the ranges before merging
    :return: List of merged ranges
    """
    if not ranges:
        raise IntRangeError('No ranges to merge')
    if len(ranges) == 1:
        yield ranges[0]
        return
    current_range = (ranges := ranges if skip_sort else sorted(ranges, key=lambda x: x[0]))[0]  # Start with the first range
    for start, end in ranges[1:]:  # Iterate through the ranges
        if start - tolerance <= current_range[1]:  # Overlap, merge the ranges
            current_range = (current_range[0], max(current_range[1], end))
        else:  # No overlap, add the current range to the merged list and start a new range
            yield current_range  # Yield the current range
            current_range = (start, end)   # Start a new range
    yield current_range  # Yield the last range


def range_overlap(range1: tuple[int, int], range2: tuple[int, int], skip_sort: bool = False) -> int:
    """
    Returns the overlap between two ranges
    :param range1: Tuple of start and end positions
    :param range2: Tuple of start and end positions
    :param skip_sort: Skip sorting each range before calculating the overlap
    :return: Integer of overlap
    """
    start1, end1 = range1 if skip_sort else sorted(range1)
    start2, end2 = range2 if skip_sort else sorted(range2)
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)
