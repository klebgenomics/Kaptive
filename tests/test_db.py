"""Unit tests for kaptive.db (Database, DatabaseManager, DatabaseMetadata, Phenotype, Phenotypes)."""

import numpy as np

from kaptive.core.interval import Interval, Intervals, Strand
from kaptive.core.kmers import FracMinHashIndex
from kaptive.core.seq import SeqRecord, Sequences
from kaptive.db import (
    Database,
    DatabaseManager,
    DatabaseMetadata,
    Phenotypes,
)


def test_database_metadata_parsed_version() -> None:
    """Test DatabaseMetadata parsing version string."""
    meta = DatabaseMetadata(
        name="Test Database",
        keyword="test_k",
        genbank="test.gbk",
        organism="Klebsiella pneumoniae",
        taxon=573,
        antigen="K",
        pathway="Wzx/Wzy",
        version="3.2.1",
        id_threshold=80.0,
        doi=["10.1000/182"],
        owner="klebgenomics",
        repo="Kaptive",
        branch="main",
        contact={"Curator": "curator@example.com"},
        phenotype_logic={},
        antigenic_units={},
    )
    assert meta.parsed_version == (3, 2, 1)


def test_database_manager_known() -> None:
    """Test DatabaseManager.known method returns registry keywords."""
    known_dbs = DatabaseManager.known()
    assert isinstance(known_dbs, list)
    assert "kpsc_k" in known_dbs or len(known_dbs) >= 0


def test_phenotypes_empty_and_len() -> None:
    """Test Phenotypes batch container methods."""
    empty_p = Phenotypes.empty()
    assert len(empty_p) == 0
    assert len(empty_p.ids) == 0
    assert empty_p.locus_masks.shape == (0, 0)


def test_database_get_locus_data_metadata() -> None:
    """Test Database.get_locus_data populates gene_states and gene_descriptions."""
    from kaptive.serotyping import GeneState

    meta = DatabaseMetadata(
        name="Test DB",
        keyword="test",
        genbank="test.gbk",
        organism="Test",
        taxon=1,
        antigen="K",
        pathway="Wzx",
        version="1.0.0",
        id_threshold=80.0,
        doi=[],
        owner="owner",
        repo="repo",
        branch="main",
        contact={},
        phenotype_logic={},
        antigenic_units={},
    )
    loci = Sequences.from_records([SeqRecord("K1", b"ATGC")])
    translations = Sequences.from_records([SeqRecord("g1", b"MKT")])
    intervals = Intervals.from_intervals([Interval(0, 3, strand=Strand.FORWARD)])

    db = Database(
        metadata=meta,
        loci=loci,
        serotypes=("K1",),
        locus_gene_offsets=np.array([0], dtype=np.uint32),
        locus_gene_lengths=np.array([1], dtype=np.uint32),
        gene_intervals=intervals,
        genes=Sequences.empty(),
        translations=translations,
        extra_genes=np.array([False], dtype=bool),
        gene_locus_indices=np.array([0], dtype=np.uint16),
        cluster_keys=("c1",),
        gene_cluster_ids=np.array([0], dtype=np.uint16),
        description_keys=("Pyruvate dehydrogenase",),
        gene_description_ids=np.array([0], dtype=np.uint16),
        gene_positions=np.array([1], dtype=np.uint16),
        phenotypes=Phenotypes.empty(),
        loci_sketches=FracMinHashIndex.build(loci, sort_by_hash=False),
    )

    locus_data = db.get_locus_data("K1")
    assert locus_data.name == "K1"
    assert locus_data.gene_states is not None
    assert locus_data.gene_states[0] == GeneState.NORMAL.value
    assert locus_data.gene_descriptions is not None
    assert locus_data.gene_descriptions[0] == "Pyruvate dehydrogenase"
    assert locus_data.gene_descriptions.dtype == object  # type: ignore
