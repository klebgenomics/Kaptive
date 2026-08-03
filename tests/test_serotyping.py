"""Unit tests for kaptive.serotyping (LocusPieces, SerotypingResult, GeneState)."""

import pytest
import numpy as np

from kaptive.serotyping import GeneState, LocusPieces


def test_locus_pieces_creation():
    """Test LocusPieces creation."""
    pieces = LocusPieces(
        ctg_indices=np.array([0, 1], dtype=np.uint32),
        starts=np.array([10, 100], dtype=np.int32),
        ends=np.array([50, 200], dtype=np.int32),
        strands=np.array([1, -1], dtype=np.int8),
    )

    assert len(pieces) == 2
    assert pieces.ctg_indices[0] == 0


def test_gene_state_enum():
    """Test GeneState enum values."""
    assert GeneState.NORMAL == 0
    assert GeneState.PARTIAL == 1
    assert GeneState.TRUNCATED == 2
    assert GeneState.NOVEL == 3


def test_serotyping_result_to_locus_data_metadata():
    """Test SerotypingResult.to_locus_data metadata fields."""
    from kaptive.core.seq import SeqRecord, Sequences
    from kaptive.serotyping import GeneHits, SerotypingResult

    gene_hits = GeneHits(
        gene_indices=np.array([0], dtype=np.int32),
        q_starts=np.array([0], dtype=np.int32),
        q_ends=np.array([100], dtype=np.int32),
        t_indices=np.array([0], dtype=np.uint32),
        t_starts=np.array([10], dtype=np.int32),
        t_ends=np.array([110], dtype=np.int32),
        strands=np.array([1], dtype=np.int8),
        is_expected=np.array([True], dtype=bool),
        is_inside=np.array([True], dtype=bool),
        is_extra=np.array([False], dtype=bool),
        expected_positions=np.array([1], dtype=np.int32),
        expected_strands=np.array([1], dtype=np.int8),
        gene_ids=np.array([b"g1"], dtype="S32"),
        cluster_names=np.array([b"c1"], dtype="S10"),
        product_descriptions=np.array([b"UDP-glucose dehydrogenase"], dtype="S64"),
        coverages=np.array([1.0], dtype=np.float32),
    )

    res = SerotypingResult(
        kaptive_version="3.0.0",
        database_name="test_db",
        database_version="1.0.0",
        database_organism="Test",
        database_taxon=1234,
        genome="Sample_1",
        best_locus_idx=0,
        best_locus_name="K1",
        best_locus_score=100.0,
        best_locus_completeness=1.0,
        length_discrepancy=0.0,
        locus_pieces=LocusPieces.empty(),
        gene_hits=gene_hits,
        gene_states=np.array([GeneState.TRUNCATED.value], dtype=np.int8),
        percent_identity=100.0,
        percent_coverage=100.0,
        phenotype="K1",
        typeable=True,
        missing_expected_genes=(),
        locus_seqs=Sequences.empty(),
        gene_seqs=Sequences.empty(),
        translations=Sequences.from_records([SeqRecord("g1", b"MKT")]),
        protein_identities=np.array([100.0], dtype=np.float32),
    )

    locus_data = res.to_locus_data()
    assert locus_data.name == "Sample_1"
    assert locus_data.gene_states is not None
    assert locus_data.gene_states[0] == GeneState.TRUNCATED.value
    assert locus_data.gene_descriptions is not None
    assert locus_data.gene_descriptions[0] == "UDP-glucose dehydrogenase"
    assert locus_data.gene_descriptions.dtype == object


def test_gene_hits_utf8_encoding():
    """Test GeneHits properly encodes non-ASCII UTF-8 characters like α, β, é."""
    from kaptive.serotyping import GeneHits

    hits = GeneHits(
        gene_indices=np.array([0], dtype=np.int32),
        q_starts=np.array([0], dtype=np.int32),
        q_ends=np.array([100], dtype=np.int32),
        t_indices=np.array([0], dtype=np.uint32),
        t_starts=np.array([10], dtype=np.int32),
        t_ends=np.array([110], dtype=np.int32),
        strands=np.array([1], dtype=np.int8),
        is_expected=np.array([True], dtype=bool),
        is_inside=np.array([True], dtype=bool),
        is_extra=np.array([False], dtype=bool),
        expected_positions=np.array([1], dtype=np.int32),
        expected_strands=np.array([1], dtype=np.int8),
        gene_ids=["gene_α"],
        cluster_names=["cluster_β"],
        product_descriptions=["synthase_é"],
        coverages=np.array([1.0], dtype=np.float32),
    )

    d = hits.to_dict()
    assert d["gene_ids"] == ["gene_α"]
    assert d["cluster_names"] == ["cluster_β"]
    assert d["product_descriptions"] == ["synthase_é"]

    reconstructed = GeneHits.from_dict(d)
    assert reconstructed.to_dict() == d


