"""
Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/kaptive-mapper

This file is part of kaptive-mapper. kaptive-mapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. kaptive-mapper is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with kaptive-mapper.
If not, see <https://www.gnu.org/licenses/>.
"""


class IntRange(object):
    """
    # TODO: KAPTIVE INTRANGE - may differ from Ryan's Verticall class - please compare.
    This class contains one or more integer ranges. Overlapping ranges will be merged together.
    It stores its ranges in a Python-like fashion where the last value in each range is
    exclusive.
    """

    def __init__(self, ranges=None):
        if not ranges:
            ranges = []
        self.ranges = []
        self.add_ranges(ranges)
        self.simplify()

    def __repr__(self):
        return str(self.ranges)

    def add_range(self, start, end):
        """Adds a single range."""
        self.add_ranges([(start, end)])

    def add_ranges(self, ranges):
        """Adds multiple ranges (list of tuples)."""
        self.ranges += ranges
        self.simplify()

    def merge_in_range(self, other):
        """Merges the other IntRange object into this one."""
        self.add_ranges(other.ranges)

    def get_total_length(self):
        """Returns the number of integers in the ranges."""
        return sum([x[1] - x[0] for x in self.ranges])

    def simplify(self):
        """Collapses overlapping ranges together."""
        fixed_ranges = []
        for int_range in self.ranges:
            if int_range[0] > int_range[1]:
                fixed_ranges.append((int_range[1], int_range[0]))
            elif int_range[0] < int_range[1]:
                fixed_ranges.append(int_range)
        starts_ends = [(x[0], 1) for x in fixed_ranges]
        starts_ends += [(x[1], -1) for x in fixed_ranges]
        starts_ends.sort(key=lambda z: z[0])
        current_sum = 0
        cumulative_sum = []
        for start_end in starts_ends:
            current_sum += start_end[1]
            cumulative_sum.append((start_end[0], current_sum))
        prev_depth = 0
        start = 0
        combined = []
        for pos, depth in cumulative_sum:
            if prev_depth == 0:
                start = pos
            elif depth == 0:
                combined.append((start, pos))
            prev_depth = depth
        self.ranges = combined

    def contains(self, other):
        """Returns True if the other IntRange is entirely contained within this IntRange."""
        for other_range in other.ranges:
            other_start, other_end = other_range
            contained = False
            for this_range in self.ranges:
                this_start, this_end = this_range
                if other_start >= this_start and other_end <= this_end:
                    contained = True
                    break
            if not contained:
                return False
        return True
