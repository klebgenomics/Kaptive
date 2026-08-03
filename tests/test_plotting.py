r"""Unit tests for kaptive.plotting (GeneStyleManager, GeneGlyphPlotter, LocusBackbonePlotter, SerotypingResultPlotter, LocusComparisonPlotter)."""

import numpy as np
import plotly.graph_objects as go
import pytest

from kaptive.compare import LocusComparisonEdges, LocusComparisons
from kaptive.core.interval import Intervals
from kaptive.plotting import (
    BasePlotter,
    GeneGlyphPlotter,
    GeneStyleManager,
    LocusBackbonePlotter,
    LocusComparisonPlotter,
    SerotypingResultPlotter,
)
from kaptive.serotyping import GeneState


def test_gene_style_manager():
    """Test GeneStyleManager state mapping and line styling."""
    assert GeneStyleManager.STATE_NAMES[GeneState.NORMAL.value] == "NORMAL"
    assert GeneStyleManager.STATE_NAMES[GeneState.PARTIAL.value] == "PARTIAL"
    assert GeneStyleManager.STATE_NAMES[GeneState.TRUNCATED.value] == "TRUNCATED"
    assert GeneStyleManager.STATE_NAMES[GeneState.NOVEL.value] == "NOVEL"

    normal_style = GeneStyleManager.get_state_style(GeneState.NORMAL, dark_mode=False)
    assert normal_style["opacity"] == 1.0
    assert normal_style["line"]["dash"] == "solid"
    assert normal_style["line"]["color"] == "black"

    truncated_style = GeneStyleManager.get_state_style(GeneState.TRUNCATED, dark_mode=False)
    assert truncated_style["opacity"] == 0.7
    assert truncated_style["line"]["dash"] == "dash"
    assert truncated_style["line"]["color"] == "orange"

    partial_style = GeneStyleManager.get_state_style(GeneState.PARTIAL, dark_mode=False)
    assert partial_style["opacity"] == 0.7
    assert partial_style["line"]["dash"] == "dash"
    assert partial_style["line"]["color"] == "red"

    novel_style = GeneStyleManager.get_state_style(GeneState.NOVEL, dark_mode=True)
    assert novel_style["opacity"] == 0.5
    assert novel_style["line"]["dash"] == "dot"
    assert novel_style["line"]["color"] == "white"


def test_gene_glyph_plotter_coordinate_generation():
    """Test GeneGlyphPlotter vectorized coordinate generation for arrows and rectangles."""
    starts = np.array([100, 500], dtype=np.float64)
    ends = np.array([400, 800], dtype=np.float64)
    strands = np.array([1, -1], dtype=np.int8)
    y_positions = np.array([0.0, 0.0], dtype=np.float64)
    head_lens = np.array([90.0, 90.0], dtype=np.float64)

    # Arrow glyphs
    x_coords, y_coords = GeneGlyphPlotter._generate_gene_coordinates(
        starts, ends, strands, y_positions, head_lens, is_rect=None
    )
    assert x_coords.shape == (2, 9)
    assert y_coords.shape == (2, 9)

    # Rectangle glyphs for truncated/partial
    is_rect = np.array([True, False], dtype=bool)
    x_rect, y_rect = GeneGlyphPlotter._generate_gene_coordinates(
        starts, ends, strands, y_positions, head_lens, is_rect=is_rect
    )
    assert x_rect.shape == (2, 9)
    assert not np.isnan(x_rect[0, 0])
    assert np.isnan(x_rect[0, 5])  # Rectangles use 5 vertices


def test_locus_backbone_plotter():
    """Test LocusBackbonePlotter trace rendering onto figure."""
    fig = go.Figure()
    LocusBackbonePlotter.render(
        fig=fig,
        x_coords=[0, 1000, None, 1200, 2000],
        y_coords=[0, 0, None, 0, 0],
        line_width=3,
        dark_mode=False,
    )
    assert len(fig.data) == 1
    assert fig.data[0].mode == "lines"
    assert fig.data[0].line.width == 3
    assert fig.data[0].line.color == "gray"


def test_locus_comparison_plotter_gene_states_and_tooltips():
    """Verify LocusComparisonPlotter renders gene states and hover tooltips accurately."""
    edges = LocusComparisonEdges.empty()
    locus_names = ("Locus_A", "Locus_B")
    locus_lengths = np.array([2, 2], dtype=np.int32)
    locus_offsets = np.array([0, 2], dtype=np.int32)
    gene_names = np.array(["geneA", "geneB", "geneC", "geneD"], dtype=object)
    gene_descriptions = np.array(
        ["Capsule synthesis A", "Capsule synthesis B", "Transporter C", "Polymerase D"],
        dtype=object,
    )
    gene_intervals = Intervals(
        starts=np.array([100, 500, 200, 700], dtype=np.int32),
        ends=np.array([400, 800, 600, 1000], dtype=np.int32),
        strands=np.array([1, -1, 1, 1], dtype=np.int8),
    )
    gene_states = np.array(
        [
            GeneState.NORMAL.value,
            GeneState.TRUNCATED.value,
            GeneState.PARTIAL.value,
            GeneState.NOVEL.value,
        ],
        dtype=np.int8,
    )

    comparisons = LocusComparisons(
        edges=edges,
        locus_names=locus_names,
        locus_lengths=locus_lengths,
        locus_offsets=locus_offsets,
        gene_names=gene_names,
        gene_descriptions=gene_descriptions,
        gene_intervals=gene_intervals,
        gene_states=gene_states,
    )

    plotter = LocusComparisonPlotter()
    fig = plotter(comparisons, align_loci=False)

    assert isinstance(fig, go.Figure)
    assert len(fig.data) > 0

    # Verify line dash styles present in figure traces (for NORMAL, TRUNCATED, PARTIAL, NOVEL)
    dash_styles = [
        t.line.dash for t in fig.data if hasattr(t, "line") and t.line and hasattr(t.line, "dash") and t.line.dash
    ]
    assert "solid" in dash_styles
    assert "dash" in dash_styles
    assert "dot" in dash_styles

    # Collect hover tooltip texts
    hover_texts = []
    for t in fig.data:
        if hasattr(t, "text") and t.text:
            if isinstance(t.text, (list, tuple, np.ndarray)):
                hover_texts.extend(t.text)

    # Verify Product: and State: in hover tooltips
    assert any("<b>Product:</b> Capsule synthesis A" in txt and "<b>State:</b> NORMAL" in txt for txt in hover_texts)
    assert any("<b>Product:</b> Capsule synthesis B" in txt and "<b>State:</b> TRUNCATED" in txt for txt in hover_texts)
    assert any("<b>Product:</b> Transporter C" in txt and "<b>State:</b> PARTIAL" in txt for txt in hover_texts)
    assert any("<b>Product:</b> Polymerase D" in txt and "<b>State:</b> NOVEL" in txt for txt in hover_texts)


def test_serotyping_result_plotter_empty():
    """Verify SerotypingResultPlotter handles empty gene hits gracefully."""
    plotter = SerotypingResultPlotter()
    # BasePlotter _get_cluster_color and layout tests
    assert BasePlotter._get_cluster_color("test_cluster").startswith("#")
    assert BasePlotter._hex_to_rgb("#FF0000") == (255, 0, 0)
