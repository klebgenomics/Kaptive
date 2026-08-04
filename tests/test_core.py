"""Unit tests for kaptive.core (interval, seq, kmers, pairwise, collections, genome)."""

import numpy as np

from kaptive.core.genome import GenomeAssembly
from kaptive.core.interval import Interval, Intervals, Strand
from kaptive.core.kmers import FracMinHashIndex, RandstrobeIndex
from kaptive.core.seq import SeqRecord, Sequences


def test_strand_coercion() -> None:
    """Test Strand enum creation and string coercion."""
    assert Strand("+") == Strand.FORWARD
    assert Strand("-") == Strand.REVERSE
    assert Strand(".") == Strand.UNSTRANDED
    assert int(Strand.FORWARD) == 1
    assert int(Strand.REVERSE) == -1
    assert int(Strand.UNSTRANDED) == 0


def test_interval_and_intervals() -> None:
    """Test Interval object and vectorized Intervals collection."""
    iv1 = Interval(10, 50, strand=Strand.FORWARD)
    iv2 = Interval(30, 80, strand=Strand.FORWARD)

    assert len(iv1) == 40
    assert iv1.start == 10
    assert iv1.end == 50
    assert 35 in iv1
    assert (5 in iv1) is False
    assert 5 not in iv1

    merged = iv1 + iv2
    assert merged.start == 10
    assert merged.end == 80

    shifted = iv1.shift(5)
    assert shifted.start == 15
    assert shifted.end == 55

    ivs = Intervals.from_intervals([iv1, iv2])
    assert len(ivs) == 2
    assert np.array_equal(ivs.starts, np.array([10, 30]))
    assert np.array_equal(ivs.ends, np.array([50, 80]))

    sliced = ivs[1]
    assert sliced.start == 30

    empty_ivs = Intervals.empty()
    assert len(empty_ivs) == 0

    concatenated = Intervals.concat([ivs, empty_ivs])
    assert len(concatenated) == 2


def test_sequences_and_seqrecord() -> None:
    """Test SeqRecord and vectorized Sequences collection."""
    rec1 = SeqRecord(id="seq1", seq=b"ATGCGTACTA")
    rec2 = SeqRecord(id="seq2", seq=b"ATGAAATAG")

    assert len(rec1) == 10
    assert rec1.id == "seq1"

    seqs = Sequences.from_records([rec1, rec2])
    assert len(seqs) == 2
    assert seqs.ids == ("seq1", "seq2")

    extracted = seqs[0]
    assert extracted.id == "seq1"
    assert extracted.seq == b"ATGCGTACTA"

    empty_seqs = Sequences.empty()
    assert len(empty_seqs) == 0

    concatenated = Sequences.concat([seqs, empty_seqs])
    assert len(concatenated) == 2


def test_kmers_indexing() -> None:
    """Test FracMinHashIndex and RandstrobeIndex k-mer construction."""
    seqs = Sequences.from_records(
        [
            SeqRecord(id="s1", seq=b"ATGCGTACTACGATCGATCGATCG"),
            SeqRecord(id="s2", seq=b"CGATCGATCGATCGATCGATCGAT"),
        ]
    )

    minhash = FracMinHashIndex.build(seqs, k=5, scaled=2)
    assert minhash.n_seqs == 2
    assert len(minhash) > 0

    strobemer = RandstrobeIndex.build(seqs)
    assert strobemer.n_seqs == 2


def test_genome_assembly() -> None:
    """Test GenomeAssembly container properties."""
    seqs = Sequences.from_records([SeqRecord(id="chr1", seq=b"ATGCGT")])
    assembly = GenomeAssembly(id="asm1", contigs=seqs)

    assert assembly.id == "asm1"
    assert len(assembly) == 6
    assert assembly.id_map == {"chr1": 0}
