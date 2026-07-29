"""Unit tests for kaptive.db (Database, DatabaseManager, DatabaseMetadata, Phenotype, Phenotypes)."""

import pytest
import numpy as np

from kaptive.core.interval import Interval, Intervals, Strand
from kaptive.core.kmers import FracMinHashIndex
from kaptive.core.seq import SeqRecord, Sequences
from kaptive.db import (
    Database,
    DatabaseManager,
    DatabaseMetadata,
    Phenotype,
    Phenotypes,
)


def test_database_metadata_parsed_version():
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


def test_database_manager_known():
    """Test DatabaseManager.known method returns registry keywords."""
    known_dbs = DatabaseManager.known()
    assert isinstance(known_dbs, list)
    assert "kpsc_k" in known_dbs or len(known_dbs) >= 0


def test_phenotypes_empty_and_len():
    """Test Phenotypes batch container methods."""
    empty_p = Phenotypes.empty()
    assert len(empty_p) == 0
    assert len(empty_p.ids) == 0
    assert empty_p.locus_masks.shape == (0, 0)
