"""Unit tests for kaptive.compare (LocusData, LocusComparisons, LocusComparisonEdges, LocusComparator)."""

import pytest
import numpy as np

from kaptive.compare import LocusComparisonEdges, LocusComparisons, LocusData
from kaptive.core.interval import Interval, Intervals, Strand
from kaptive.core.seq import SeqRecord, Sequences


def test_locus_data_creation():
    """Test LocusData creation and property evaluation."""
    proteins = Sequences.from_records(
        [
            SeqRecord(id="p1", seq=b"MKTLLILAV"),
            SeqRecord(id="p2", seq=b"MSKGEELFTG"),
        ]
    )
    intervals = Intervals.from_intervals(
        [
            Interval(10, 100, strand=Strand.FORWARD),
            Interval(200, 500, strand=Strand.REVERSE),
        ]
    )

    locus_data = LocusData(
        name="Locus_1",
        proteins=proteins,
        backbone=intervals,
    )

    assert locus_data.name == "Locus_1"
    assert len(locus_data.proteins) == 2
    assert len(locus_data.backbone) == 2
