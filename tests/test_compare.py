"""Unit tests for kaptive.compare (LocusData, LocusComparisons, LocusComparisonEdges, LocusComparator)."""

import numpy as np
import pytest

from kaptive.compare import LocusData
from kaptive.core.interval import Interval, Intervals, Strand
from kaptive.core.seq import SeqRecord, Sequences


def test_locus_data_creation() -> None:
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
    assert locus_data.gene_states is None
    assert locus_data.gene_descriptions is None


def test_locus_comparator_metadata_propagation() -> None:
    """Test LocusComparator concatenating gene_states and gene_descriptions."""
    from kaptive.compare import LocusComparator
    from kaptive.serotyping import GeneState

    p1 = Sequences.from_records([SeqRecord(id="g1", seq=b"MKTLLILAV")])
    b1 = Intervals.from_intervals([Interval(10, 100, strand=Strand.FORWARD)])
    l1 = LocusData(
        name="Locus_1",
        proteins=p1,
        backbone=b1,
        gene_states=np.array([GeneState.TRUNCATED.value], dtype=np.int8),
        gene_descriptions=np.array(["glycosyltransferase"], dtype=object),
    )

    p2 = Sequences.from_records([SeqRecord(id="g2", seq=b"MSKGEELFTG")])
    b2 = Intervals.from_intervals([Interval(200, 500, strand=Strand.FORWARD)])
    l2 = LocusData(
        name="Locus_2",
        proteins=p2,
        backbone=b2,
        # gene_states and gene_descriptions left None to test defaulting
    )

    comparator = LocusComparator()
    comp = comparator([l1, l2])

    assert len(comp.gene_states) == 2
    assert comp.gene_states[0] == GeneState.TRUNCATED.value
    assert comp.gene_states[1] == GeneState.NORMAL.value
    assert comp.gene_descriptions[0] == "glycosyltransferase"
    assert comp.gene_descriptions[1] == ""


def test_locus_comparator_bytes_decoding() -> None:
    """Test LocusComparator properly decodes bytes and object arrays containing bytes to UTF-8 strings."""
    from kaptive.compare import LocusComparator

    p1 = Sequences.from_records(
        [
            SeqRecord(id="g1", seq=b"MKTLLILAV"),
            SeqRecord(id="g2", seq=b"MSKGEELFTG"),
        ]
    )
    b1 = Intervals.from_intervals(
        [
            Interval(10, 100, strand=Strand.FORWARD),
            Interval(200, 500, strand=Strand.FORWARD),
        ]
    )

    # Case 1: Bytes array (S dtype)
    l1 = LocusData(
        name="Locus_1",
        proteins=p1,
        backbone=b1,
        gene_descriptions=np.array([b"wzi_\xce\xb1", b"wza_\xce\xb2"]),
    )

    # Case 2: Object array containing bytes
    l2 = LocusData(
        name="Locus_2",
        proteins=p1,
        backbone=b1,
        gene_descriptions=np.array([b"wzi_\xce\xb1", b"wza_\xce\xb2"], dtype=object),
    )

    comparator = LocusComparator()
    comp = comparator([l1, l2])

    assert comp.gene_descriptions.dtype == object
    assert list(comp.gene_descriptions) == ["wzi_α", "wza_β", "wzi_α", "wza_β"]


def test_locus_comparator_shape_validation() -> None:
    """Test LocusComparator raises ValueError when array shapes mismatch len(proteins)."""
    from kaptive.compare import LocusComparator

    p1 = Sequences.from_records([SeqRecord(id="g1", seq=b"MKTLLILAV")])
    b1 = Intervals.from_intervals([Interval(10, 100, strand=Strand.FORWARD)])

    # Mismatched gene_descriptions length
    l_bad_desc = LocusData(
        name="Locus_BadDesc",
        proteins=p1,
        backbone=b1,
        gene_descriptions=["desc1", "desc2"],
    )

    comparator = LocusComparator()
    with pytest.raises(ValueError, match="gene_descriptions length"):
        comparator([l_bad_desc])

    # Mismatched gene_states length
    l_bad_state = LocusData(
        name="Locus_BadState",
        proteins=p1,
        backbone=b1,
        gene_states=np.array([0, 1], dtype=np.int8),
    )

    with pytest.raises(ValueError, match="gene_states length"):
        comparator([l_bad_state])

    # Mismatched backbone length
    b_bad = Intervals.from_intervals(
        [
            Interval(10, 100, strand=Strand.FORWARD),
            Interval(200, 500, strand=Strand.FORWARD),
        ]
    )
    l_bad_bb = LocusData(
        name="Locus_BadBB",
        proteins=p1,
        backbone=b_bad,
    )

    with pytest.raises(ValueError, match="backbone length"):
        comparator([l_bad_bb])
