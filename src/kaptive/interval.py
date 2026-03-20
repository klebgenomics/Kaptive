from typing import Generator
from operator import itemgetter

# Functions ------------------------------------------------------------------------------------------------------------
def merge_ranges(ranges: list[tuple[int | float, int | float]], tolerance: int | float = 0, skip_sort: bool = False
                 ) -> Generator[tuple[int | float, int | float], None, None]:
    """
    Merge overlapping ranges
    :param ranges: List of tuples of start and end positions
    :param tolerance: Integer or float of tolerance for merging ranges
    :param skip_sort: Skip sorting the ranges before merging
    :return: List of merged ranges
    """
    if not ranges:
        return None
    if len(ranges) == 1:
        yield ranges[0]
        return None
    current_range = (ranges := ranges if skip_sort else sorted(ranges, key=itemgetter(0)))[0]  # Start with the first range
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

