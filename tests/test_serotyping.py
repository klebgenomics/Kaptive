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
