r"""Interactive Plotly visualization components for Kaptive serotyping and comparison maps.

Provides modular OOP plotting components for:
- Single-locus serotyping results ([`SerotypingResultPlotter`][kaptive.plotting.SerotypingResultPlotter])
- Multi-locus comparative synteny maps ([`LocusComparisonPlotter`][kaptive.plotting.LocusComparisonPlotter])

Classes:
    [`BasePlotter`][kaptive.plotting.BasePlotter]: Abstract base class providing color palette and layout helpers.
    [`GeneStyleManager`][kaptive.plotting.GeneStyleManager]: Manager for gene state visual styles and state names.
    [`GeneGlyphPlotter`][kaptive.plotting.GeneGlyphPlotter]: Component for rendering gene polygon glyphs and hover markers.
    [`LocusBackbonePlotter`][kaptive.plotting.LocusBackbonePlotter]: Component for rendering locus backbone lines.
    [`SerotypingResultPlotter`][kaptive.plotting.SerotypingResultPlotter]: Plotter for single-locus serotyping results.
    [`LocusComparisonPlotter`][kaptive.plotting.LocusComparisonPlotter]: Plotter for multi-locus comparison graphs.
"""

import zlib
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Sequence
from typing import Any, ClassVar

import numba
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from kaptive.compare import LocusComparisons
from kaptive.core.interval import Intervals
from kaptive.serotyping import GeneState, SerotypingResult


# Classes --------------------------------------------------------------------------------------------------------------
class BasePlotter(ABC):
    r"""Abstract base plotter providing common color palette and layout helpers.

    Attributes:
        PALETTE: Categorical color palette used for cluster and gene styling.
    """

    PALETTE = px.colors.qualitative.Alphabet

    @classmethod
    def _get_cluster_color(cls, cluster_id_or_name: int | str) -> str:
        r"""Get a deterministic hex color for a given cluster ID or name.

        Args:
            cluster_id_or_name (int | str): Cluster integer ID or string identifier.

        Returns:
            str: Hexadecimal color string.
        """
        if isinstance(cluster_id_or_name, (int, np.integer)):
            if cluster_id_or_name < 0:
                return "#D3D3D3"
            idx = cluster_id_or_name % len(cls.PALETTE)
        else:
            idx = zlib.crc32(str(cluster_id_or_name).encode()) % len(cls.PALETTE)
        return cls.PALETTE[idx]

    @classmethod
    def _hex_to_rgb(cls, hex_col: str) -> tuple[int, int, int]:
        r"""Convert a hex color string into an RGB integer tuple.

        Args:
            hex_col (str): Hexadecimal color string (e.g. '#FF0000').

        Returns:
            tuple[int, int, int]: RGB color tuple.
        """
        hex_col = hex_col.lstrip("#")
        if len(hex_col) == 6:
            return tuple(int(hex_col[i : i + 2], 16) for i in (0, 2, 4))
        return (128, 128, 128)

    @classmethod
    def _generate_gene_coordinates(
        cls,
        starts: np.ndarray,
        ends: np.ndarray,
        strands: np.ndarray,
        y_positions: np.ndarray,
        head_lens: np.ndarray,
        is_rect: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Generate 2D polygon vertex coordinates for drawing gene glyphs.

        Delegates to [`GeneGlyphPlotter._generate_gene_coordinates`][kaptive.plotting.GeneGlyphPlotter._generate_gene_coordinates].

        Args:
            starts (np.ndarray): Gene start coordinates.
            ends (np.ndarray): Gene end coordinates.
            strands (np.ndarray): Gene strand orientations (+1 for forward, -1 for reverse).
            y_positions (np.ndarray): Baseline Y-axis positions for genes.
            head_lens (np.ndarray): Arrowhead lengths.
            is_rect (np.ndarray | None): Boolean mask specifying rectangular rendering per gene.

        Returns:
            tuple[np.ndarray, np.ndarray]: X and Y coordinate matrices of shape (N, 9).
        """
        return GeneGlyphPlotter._generate_gene_coordinates(
            starts, ends, strands, y_positions, head_lens, is_rect=is_rect
        )

    @classmethod
    def _apply_standard_layout(
        cls,
        fig: go.Figure,
        title: str | None = None,
        xaxis: dict[str, Any] | None = None,
        yaxis: dict[str, Any] | None = None,
        dark_mode: bool = False,
    ) -> None:
        r"""Apply standard Plotly figure layout updates.

        Args:
            fig (go.Figure): Target Plotly figure.
            title (str | None): Optional figure title.
            xaxis (dict[str, Any] | None): Optional X-axis configuration dictionary.
            yaxis (dict[str, Any] | None): Optional Y-axis configuration dictionary.
            dark_mode (bool): Whether dark mode styling should be applied.
        """
        layout_kwargs: dict[str, Any] = dict(
            template="plotly_dark" if dark_mode else "plotly_white",
            paper_bgcolor="rgba(0,0,0,0)" if dark_mode else "white",
            plot_bgcolor="rgba(0,0,0,0)" if dark_mode else "white",
            hovermode="closest",
        )
        if title is not None:
            layout_kwargs["title"] = title
        if xaxis is not None:
            layout_kwargs["xaxis"] = xaxis
        if yaxis is not None:
            layout_kwargs["yaxis"] = yaxis

        fig.update_layout(**layout_kwargs)

    @abstractmethod
    def __call__(self, *args: Any, **kwargs: Any) -> go.Figure:
        r"""Abstract call interface for rendering a Plotly figure.

        Returns:
            go.Figure: Interactive Plotly figure.
        """
        pass


class GeneStyleManager:
    r"""Manages visual attributes (colors, opacities, stroke styles) for gene states.

    Attributes:
        STATE_NAMES (ClassVar[dict[int, str]]): Mapping from GeneState integer value to human-readable string.
    """

    STATE_NAMES: ClassVar[dict[int, str]] = {s.value: s.name for s in GeneState}

    @classmethod
    def get_state_style(cls, state: int | GeneState | str, dark_mode: bool = False) -> dict[str, Any]:
        r"""Get Plotly line and opacity style dict for a given gene state.

        Args:
            state (int | GeneState | str): GeneState value, integer, or name.
            dark_mode (bool): Whether dark theme line colors should be used.

        Returns:
            dict[str, Any]: Dictionary containing 'opacity' (float) and 'line' (dict).
        """
        if isinstance(state, GeneState):
            st_val = state.value
        elif isinstance(state, str):
            try:
                st_val = GeneState[state].value
            except KeyError:
                st_val = int(state) if state.isdigit() else 0
        else:
            st_val = int(state)

        line_color = "white" if dark_mode else "black"
        styles = {
            GeneState.NORMAL.value: dict(opacity=1.0, line=dict(width=1, dash="solid", color=line_color)),
            GeneState.NOVEL.value: dict(opacity=0.5, line=dict(width=1, dash="dot", color=line_color)),
            GeneState.PARTIAL.value: dict(opacity=0.7, line=dict(width=1, dash="dash", color="red")),
            GeneState.TRUNCATED.value: dict(opacity=0.7, line=dict(width=1, dash="dash", color="orange")),
        }
        return styles.get(st_val, dict(opacity=1.0, line=dict(width=1, dash="solid", color=line_color)))


class GeneGlyphPlotter:
    r"""Component for generating and rendering gene polygon glyphs and hover markers in Plotly figures."""

    @classmethod
    def _generate_gene_coordinates(
        cls,
        starts: np.ndarray,
        ends: np.ndarray,
        strands: np.ndarray,
        y_positions: np.ndarray,
        head_lens: np.ndarray,
        is_rect: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Generate 2D polygon vertex coordinates for drawing gene glyphs.

        Computes 8-vertex directional arrow polygons (or 5-vertex rectangular boxes if `is_rect` is set)
        for vectorized rendering in Plotly scatter traces.

        Args:
            starts (np.ndarray): Gene start coordinates.
            ends (np.ndarray): Gene end coordinates.
            strands (np.ndarray): Gene strand orientations (+1 for forward, -1 for reverse).
            y_positions (np.ndarray): Baseline Y-axis positions for genes.
            head_lens (np.ndarray): Arrowhead lengths.
            is_rect (np.ndarray | None): Boolean mask specifying rectangular rendering per gene.

        Returns:
            tuple[np.ndarray, np.ndarray]: X and Y coordinate matrices of shape (N, 9).
        """
        N = len(starts)
        x_coords = np.full((N, 9), np.nan, dtype=np.float64)
        y_coords = np.full((N, 9), np.nan, dtype=np.float64)

        if is_rect is not None and np.any(is_rect):
            s = starts[is_rect]
            e = ends[is_rect]
            x_coords[is_rect, 0] = s
            x_coords[is_rect, 1] = e
            x_coords[is_rect, 2] = e
            x_coords[is_rect, 3] = s
            x_coords[is_rect, 4] = s
            y = y_positions[is_rect]
            y_coords[is_rect, 0] = y - 0.2
            y_coords[is_rect, 1] = y - 0.2
            y_coords[is_rect, 2] = y + 0.2
            y_coords[is_rect, 3] = y + 0.2
            y_coords[is_rect, 4] = y - 0.2

        if is_rect is not None:
            is_arrow = ~is_rect
        else:
            is_arrow = np.ones(N, dtype=bool)

        is_fwd = is_arrow & (strands >= 0)
        is_rev = is_arrow & (strands < 0)

        if np.any(is_fwd):
            s, e, hl, y = starts[is_fwd], ends[is_fwd], head_lens[is_fwd], y_positions[is_fwd]
            x_coords[is_fwd, 0] = s
            x_coords[is_fwd, 1] = np.maximum(s, e - hl)
            x_coords[is_fwd, 2] = np.maximum(s, e - hl)
            x_coords[is_fwd, 3] = e
            x_coords[is_fwd, 4] = np.maximum(s, e - hl)
            x_coords[is_fwd, 5] = np.maximum(s, e - hl)
            x_coords[is_fwd, 6] = s
            x_coords[is_fwd, 7] = s

            y_coords[is_fwd, 0] = y - 0.2
            y_coords[is_fwd, 1] = y - 0.2
            y_coords[is_fwd, 2] = y - 0.3
            y_coords[is_fwd, 3] = y + 0.0
            y_coords[is_fwd, 4] = y + 0.3
            y_coords[is_fwd, 5] = y + 0.2
            y_coords[is_fwd, 6] = y + 0.2
            y_coords[is_fwd, 7] = y - 0.2

        if np.any(is_rev):
            s, e, hl, y = starts[is_rev], ends[is_rev], head_lens[is_rev], y_positions[is_rev]
            x_coords[is_rev, 0] = e
            x_coords[is_rev, 1] = np.minimum(e, s + hl)
            x_coords[is_rev, 2] = np.minimum(e, s + hl)
            x_coords[is_rev, 3] = s
            x_coords[is_rev, 4] = np.minimum(e, s + hl)
            x_coords[is_rev, 5] = np.minimum(e, s + hl)
            x_coords[is_rev, 6] = e
            x_coords[is_rev, 7] = e

            y_coords[is_rev, 0] = y - 0.2
            y_coords[is_rev, 1] = y - 0.2
            y_coords[is_rev, 2] = y - 0.3
            y_coords[is_rev, 3] = y + 0.0
            y_coords[is_rev, 4] = y + 0.3
            y_coords[is_rev, 5] = y + 0.2
            y_coords[is_rev, 6] = y + 0.2
            y_coords[is_rev, 7] = y - 0.2

        return x_coords, y_coords

    @classmethod
    def render(
        cls,
        fig: go.Figure,
        starts: np.ndarray,
        ends: np.ndarray,
        strands: np.ndarray,
        y_positions: np.ndarray,
        cluster_keys: Sequence[Any],
        states: np.ndarray | None,
        hover_texts: Sequence[str],
        legend_names: Sequence[str] | None = None,
        dark_mode: bool = False,
        seen_clusters: set[str] | None = None,
        showlegend: bool = True,
    ) -> set[str]:
        r"""Render gene glyph polygons and invisible hover markers onto a Plotly figure.

        Args:
            fig (go.Figure): Target Plotly figure.
            starts (np.ndarray): Gene start X-coordinates.
            ends (np.ndarray): Gene end X-coordinates.
            strands (np.ndarray): Gene strands (+1 or -1).
            y_positions (np.ndarray): Baseline Y-coordinates.
            cluster_keys (Sequence[Any]): Cluster identifiers for color lookup.
            states (np.ndarray | None): Array of GeneState values per gene.
            hover_texts (Sequence[str]): Formatted hover text strings per gene.
            legend_names (Sequence[str] | None): Optional custom legend entry names per gene.
            dark_mode (bool): Whether dark theme styling is enabled.
            seen_clusters (set[str] | None): Set of cluster legend names already rendered.
            showlegend (bool): Whether legend entries should be displayed.

        Returns:
            set[str]: Updated set of cluster legend names.
        """
        N = len(starts)
        if N == 0:
            return seen_clusters or set()

        if seen_clusters is None:
            seen_clusters = set()

        if states is None:
            states = np.zeros(N, dtype=np.int8)

        # Shape logic: rectangle for PARTIAL/TRUNCATED, arrow for NORMAL/NOVEL
        is_rect = ~((states == GeneState.NORMAL.value) | (states == GeneState.NOVEL.value))

        lengths = ends - starts
        head_lens = np.minimum(lengths * 0.3, 400.0)

        x_coords, y_coords = cls._generate_gene_coordinates(
            starts, ends, strands, y_positions, head_lens, is_rect=is_rect
        )

        if legend_names is None:
            leg_list = []
            for ck in cluster_keys:
                if isinstance(ck, bytes):
                    ck_val: Any = ck.decode()
                else:
                    ck_val = ck
                if isinstance(ck_val, (int, np.integer)):
                    leg_list.append(f"Cluster {ck_val}" if ck_val >= 0 else "Unique")
                else:
                    leg_list.append(str(ck_val))
        else:
            leg_list = list(legend_names)

        # Group genes by (legend_name, cluster_key, state)
        traces_data = defaultdict(list)
        for i in range(N):
            ck = cluster_keys[i]
            if isinstance(ck, bytes):
                ck = ck.decode()
            traces_data[(leg_list[i], ck, states[i])].append(i)

        for (leg_name, ck, st), idxs in traces_data.items():
            x_vals = x_coords[idxs].flatten().tolist()
            y_vals = y_coords[idxs].flatten().tolist()

            color = BasePlotter._get_cluster_color(ck)
            style = GeneStyleManager.get_state_style(st, dark_mode=dark_mode)

            show_leg = showlegend and (leg_name not in seen_clusters)
            seen_clusters.add(leg_name)

            fig.add_trace(
                go.Scatter(
                    x=x_vals,
                    y=y_vals,
                    mode="lines",
                    fill="toself",
                    fillcolor=color,
                    line=style["line"],
                    opacity=style["opacity"],
                    name=leg_name,
                    legendgroup=leg_name,
                    showlegend=show_leg,
                    hoverinfo="none",
                )
            )

        # Add invisible hover markers
        cx = ((starts + ends) / 2.0).tolist()
        cy = y_positions.tolist()

        hover_colors = []
        for ck in cluster_keys:
            if isinstance(ck, bytes):
                ck_val_c: Any = ck.decode()
            else:
                ck_val_c = ck
            hover_colors.append(BasePlotter._get_cluster_color(ck_val_c))

        fig.add_trace(
            go.Scatter(
                x=cx,
                y=cy,
                mode="markers",
                marker=dict(size=1, color="rgba(0,0,0,0)"),
                text=list(hover_texts),
                hoverinfo="text",
                hoverlabel=dict(bgcolor=hover_colors),
                showlegend=False,
            )
        )

        return seen_clusters


class LocusBackbonePlotter:
    r"""Component for rendering horizontal locus backbone lines in Plotly figures."""

    @classmethod
    def render(
        cls,
        fig: go.Figure,
        x_coords: list[float | None] | np.ndarray,
        y_coords: list[float | None] | np.ndarray,
        line_width: int = 2,
        dark_mode: bool = False,
    ) -> None:
        r"""Render horizontal backbone line segments onto a Plotly figure.

        Args:
            fig (go.Figure): Target Plotly figure.
            x_coords (list[float | None] | np.ndarray): X coordinates of backbone segments.
            y_coords (list[float | None] | np.ndarray): Y coordinates of backbone segments.
            line_width (int): Stroke width of the backbone line. Defaults to 2.
            dark_mode (bool): Whether dark mode theme colors should be used. Defaults to False.
        """
        line_color = "white" if dark_mode else "gray"
        fig.add_trace(
            go.Scatter(
                x=list(x_coords),
                y=list(y_coords),
                mode="lines",
                line=dict(color=line_color, width=line_width),
                hoverinfo="none",
                showlegend=False,
            )
        )


class SerotypingResultPlotter(BasePlotter):
    r"""Plotter for visualising single-locus serotyping results.

    Generates interactive Plotly diagrams illustrating gene hit locations, locus contig
    backbones, gene state annotations, and detailed hover metadata.
    """

    STATE_NAMES: ClassVar[dict[int, str]] = GeneStyleManager.STATE_NAMES

    @classmethod
    def _get_state_styles(cls, dark_mode: bool = False) -> dict[int, dict[str, Any]]:
        r"""Get line and opacity style mappings for gene states.

        Args:
            dark_mode (bool): Whether dark mode palette styles should be used. Defaults to False.

        Returns:
            dict[int, dict[str, Any]]: Dictionary mapping [`GeneState`][kaptive.serotyping.GeneState] values
                to Plotly line and opacity style dicts.
        """
        return {st: GeneStyleManager.get_state_style(st, dark_mode) for st in cls.STATE_NAMES}

    def __call__(self, result: SerotypingResult, dark_mode: bool = False) -> go.Figure:
        r"""Render an interactive Plotly diagram for a serotyping result.

        Args:
            result (SerotypingResult): The serotyping result to visualize.
            dark_mode (bool): Whether to render using dark theme styling. Defaults to False.

        Returns:
            go.Figure: Interactive Plotly figure displaying locus gene layout and metadata.
        """
        # Filter for genes inside the locus
        mask = result.gene_hits.is_inside
        genes = result.gene_hits[mask]
        states = result.gene_states[mask]
        identities = result.protein_identities[mask]
        expected_positions = genes.expected_positions

        N = len(genes)
        if N == 0:
            return go.Figure()

        # Determine which piece each gene belongs to
        piece_indices = np.zeros(N, dtype=np.int32)
        max_overlaps = np.zeros(N, dtype=np.int32) - 1
        n_pieces = len(result.locus_pieces)
        for i in range(n_pieces):
            p_ctg = result.locus_pieces.ctg_indices[i]
            p_s = result.locus_pieces.starts[i]
            p_e = result.locus_pieces.ends[i]

            same_ctg = genes.t_indices == p_ctg
            overlap_starts = np.maximum(genes.t_starts, p_s)
            overlap_ends = np.minimum(genes.t_ends, p_e)
            overlaps = np.maximum(0, overlap_ends - overlap_starts)

            better_overlap = same_ctg & (overlaps > max_overlaps)
            piece_indices[better_overlap] = i
            max_overlaps[better_overlap] = overlaps[better_overlap]

        # Calculate piece order based on expected positions
        piece_mean_pos = np.zeros(n_pieces, dtype=np.float64)

        for i in range(n_pieces):
            in_piece = piece_indices == i
            if not np.any(in_piece):
                continue

            expected_mask = in_piece & genes.is_expected
            if np.any(expected_mask):
                piece_mean_pos[i] = np.mean(expected_positions[expected_mask])
            else:
                piece_mean_pos[i] = np.mean(expected_positions[in_piece])

        piece_order = np.argsort(piece_mean_pos).astype(np.int32)

        # Use Intervals to vectorise the layout geometry
        genes_intervals = Intervals(genes.t_starts, genes.t_ends, genes.strands)
        arranged_intervals = genes_intervals.arrange(
            piece_indices,
            piece_order,
            result.locus_pieces.starts,
            result.locus_pieces.ends,
            result.locus_pieces.strands,
        )

        plot_starts = arranged_intervals.starts
        plot_ends = arranged_intervals.ends
        plot_strands = arranged_intervals.strands

        piece_plot_starts = np.zeros(n_pieces, dtype=np.int32)
        piece_plot_ends = np.zeros(n_pieces, dtype=np.int32)
        gap = 500
        current_x = 0
        for i in piece_order:
            p_len = result.locus_pieces.ends[i] - result.locus_pieces.starts[i]
            piece_plot_starts[i] = current_x
            piece_plot_ends[i] = current_x + p_len
            current_x += p_len + gap

        fig = go.Figure()

        # Add contig backbone via LocusBackbonePlotter
        piece_x: list[float | None] = []
        piece_y: list[float | None] = []
        tick_vals: list[float] = []
        tick_text: list[str] = []

        for i in piece_order:
            p_s = piece_plot_starts[i]
            p_e = piece_plot_ends[i]
            piece_x.extend([p_s, p_e, None])
            piece_y.extend([0, 0, None])

            ctg_name = result.locus_seqs.ids[i]
            parts = ctg_name.split("_")
            orient = result.locus_pieces.strands[i]

            if len(parts) >= 4 and parts[-1] == "0":
                name = "_".join(parts[:-3])
                start = parts[-3]
                end = parts[-2]
                if orient == 1:
                    label = f"{name}:{start}-{end}"
                else:
                    label = f"{name}:{end}-{start}"
            else:
                label = ctg_name

            tick_vals.append((p_s + p_e) / 2.0)
            tick_text.append(label)

        LocusBackbonePlotter.render(fig, piece_x, piece_y, line_width=4, dark_mode=dark_mode)

        # Prepare hover texts
        cluster_names = np.char.decode(genes.cluster_names)
        gene_ids = np.char.decode(genes.gene_ids)
        product_descs = np.char.decode(genes.product_descriptions)
        strands_text = ["+" if st > 0 else "-" for st in genes.strands]

        hover_texts = [
            f"<b>Gene:</b> {g}<br>"
            f"<b>Cluster:</b> {c}<br>"
            f"<b>Location:</b> {start}-{end} ({strand})<br>"
            f"<b>Coverage:</b> {cov:.1f}%<br>"
            f"<b>Identity:</b> {ident:.1f}%<br>"
            f"<b>State:</b> {GeneStyleManager.STATE_NAMES.get(int(s), 'NORMAL')}<br>"
            f"<b>Product:</b> {d}"
            for g, c, start, end, strand, cov, ident, s, d in zip(
                gene_ids,
                cluster_names,
                genes.t_starts,
                genes.t_ends,
                strands_text,
                genes.coverages,
                identities,
                states,
                product_descs,
            )
        ]

        y_positions = np.zeros(N, dtype=np.float64)

        # Render gene glyphs & hover markers via GeneGlyphPlotter
        GeneGlyphPlotter.render(
            fig=fig,
            starts=plot_starts,
            ends=plot_ends,
            strands=plot_strands,
            y_positions=y_positions,
            cluster_keys=cluster_names,
            states=states,
            hover_texts=hover_texts,
            legend_names=cluster_names,
            dark_mode=dark_mode,
        )

        self._apply_standard_layout(
            fig=fig,
            title=f"Serotype Result: {result.best_locus_name} ({result.phenotype})",
            xaxis=dict(
                tickmode="array",
                tickvals=tick_vals,
                ticktext=tick_text,
                showgrid=False,
                zeroline=False,
            ),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False, range=[-1, 1]),
            dark_mode=dark_mode,
        )

        return fig


class LocusComparisonPlotter(BasePlotter):
    r"""Plotter for multi-locus comparative synteny maps.

    Renders stacked multi-locus gene diagrams with homology ribbons connecting homologous
    genes across loci, offset alignments, and interactive hover tooltips.
    """

    @staticmethod
    @numba.njit(cache=True, nogil=True)
    def _connected_components(n: int, edges_u: np.ndarray, edges_v: np.ndarray) -> np.ndarray:
        r"""Compute connected components of a gene homology graph using disjoint-set union-find.

        Args:
            n (int): Total number of genes (nodes).
            edges_u (np.ndarray): 1D array of query gene indices (edges source).
            edges_v (np.ndarray): 1D array of target gene indices (edges target).

        Returns:
            np.ndarray: Component cluster ID assignment array of shape (n,).
        """
        parent = np.arange(n, dtype=np.int32)

        def find(i):
            root = i
            while root != parent[root]:
                root = parent[root]
            curr = i
            while curr != root:
                nxt = parent[curr]
                parent[curr] = root
                curr = nxt
            return root

        def union(i, j):
            root_i = find(i)
            root_j = find(j)
            if root_i != root_j:
                parent[root_i] = root_j

        for i in range(len(edges_u)):
            union(edges_u[i], edges_v[i])

        for i in range(n):
            parent[i] = find(i)

        return parent

    @staticmethod
    def _calculate_alignment_offsets(comparisons: LocusComparisons) -> np.ndarray:
        r"""Calculate horizontal X-axis offsets to align syntenic gene regions across loci.

        Args:
            comparisons (LocusComparisons): Multi-locus comparison result container.

        Returns:
            np.ndarray: 1D array of integer X-axis offsets for each locus.
        """
        n_loci = len(comparisons.locus_lengths)
        offsets = np.zeros(n_loci, dtype=np.float64)
        if n_loci <= 1:
            return offsets.astype(np.int32)

        for i in range(1, n_loci):
            mask = (comparisons.edges.query_locus_indices == i - 1) & (comparisons.edges.target_locus_indices == i)
            if not np.any(mask):
                found = False
                for prev_i in range(i - 1, -1, -1):
                    mask = (comparisons.edges.query_locus_indices == prev_i) & (
                        comparisons.edges.target_locus_indices == i
                    )
                    if np.any(mask):
                        found = True
                        break

                if not found:
                    continue

            q_idx = comparisons.edges.global_query_indices[mask]
            t_idx = comparisons.edges.global_target_indices[mask]

            prev_i = comparisons.edges.query_locus_indices[mask][0]

            q_starts = comparisons.gene_intervals.starts[q_idx]
            t_starts = comparisons.gene_intervals.starts[t_idx]

            shifts = (q_starts + offsets[prev_i]) - t_starts
            offsets[i] = np.median(shifts)

        # Normalize so the minimum coordinate is 0
        min_coords = []
        for i in range(n_loci):
            l_start = comparisons.locus_offsets[i]
            l_len = comparisons.locus_lengths[i]
            if l_len > 0:
                min_coords.append(np.min(comparisons.gene_intervals.starts[l_start : l_start + l_len]) + offsets[i])

        if min_coords:
            offsets -= np.min(min_coords)

        return offsets.astype(np.int32)

    def __call__(
        self,
        comparisons: LocusComparisons,
        offsets: Sequence[int] | None = None,
        align_loci: bool = True,
        show_all_links: bool = False,
        dark_mode: bool = False,
    ) -> go.Figure:
        r"""Render an interactive Plotly diagram for multi-locus comparative analysis.

        Args:
            comparisons (LocusComparisons): Multi-locus comparison container.
            offsets (Sequence[int] | None): Custom per-locus X-axis offsets. Defaults to None.
            align_loci (bool): Whether to calculate alignment offsets automatically. Defaults to True.
            show_all_links (bool): Whether to draw homology links between all locus pairs
                or adjacent loci only. Defaults to False.
            dark_mode (bool): Whether to render using dark theme styling. Defaults to False.

        Returns:
            go.Figure: Interactive Plotly multi-locus comparison figure.
        """
        n_loci = len(comparisons.locus_lengths)
        if n_loci == 0:
            return go.Figure()

        if offsets is None:
            if align_loci:
                calculated_offsets = self._calculate_alignment_offsets(comparisons)
            else:
                calculated_offsets = np.zeros(n_loci, dtype=np.int32)
        else:
            calculated_offsets = np.asarray(offsets, dtype=np.int32)

        total_genes = len(comparisons.gene_intervals)
        if total_genes == 0:
            return go.Figure()

        edges = comparisons.edges

        # 1. Build Graph and cluster (connected components)
        if len(edges) > 0:
            labels = self._connected_components(total_genes, edges.global_query_indices, edges.global_target_indices)

            # Find which clusters are actually non-singletons (Components with >1 nodes)
            unique, counts = np.unique(labels, return_counts=True)
            singleton_labels = unique[counts == 1]
            is_singleton = np.isin(labels, singleton_labels)
            labels[is_singleton] = -1
        else:
            labels = np.full(total_genes, -1, dtype=np.int32)

        fig = go.Figure()

        # 2. Extract spatial arrays
        global_to_locus = np.repeat(np.arange(n_loci, dtype=np.int32), comparisons.locus_lengths)

        starts = comparisons.gene_intervals.starts + calculated_offsets[global_to_locus]
        ends = comparisons.gene_intervals.ends + calculated_offsets[global_to_locus]
        strands = comparisons.gene_intervals.strands

        y_positions = -np.arange(n_loci, dtype=np.float64) * 2.0

        # 3. Plot Backbones using LocusBackbonePlotter
        backbone_x: list[float | None] = []
        backbone_y: list[float | None] = []
        for i in range(n_loci):
            l_start = comparisons.locus_offsets[i]
            l_len = comparisons.locus_lengths[i]
            if l_len > 0:
                min_x = np.min(starts[l_start : l_start + l_len])
                max_x = np.max(ends[l_start : l_start + l_len])
                backbone_x.extend([min_x, max_x, None])
                backbone_y.extend([y_positions[i], y_positions[i], None])

        if backbone_x:
            LocusBackbonePlotter.render(fig, backbone_x, backbone_y, line_width=2, dark_mode=dark_mode)

        # 4. Add Homology Links (Polygons with opacity binning)
        if len(edges) > 0:
            link_mask = np.ones(len(edges), dtype=bool)
            if not show_all_links:
                link_mask = edges.target_locus_indices == edges.query_locus_indices + 1

            q_idx = edges.global_query_indices[link_mask]
            t_idx = edges.global_target_indices[link_mask]
            pidents = edges.alignments.pidents[link_mask]

            if len(q_idx) > 0:
                # Opacity binning via digitize
                min_pid, max_pid = np.min(pidents), np.max(pidents)
                if max_pid > min_pid:
                    bins = np.linspace(min_pid, max_pid, 6)  # 5 bins
                    binned_idxs = np.digitize(pidents, bins[1:], right=True)
                else:
                    binned_idxs = np.zeros(len(pidents), dtype=np.int32)

                q_starts, q_ends = starts[q_idx], ends[q_idx]
                t_starts, t_ends = starts[t_idx], ends[t_idx]

                q_y = y_positions[global_to_locus[q_idx]] - 0.2
                t_y = y_positions[global_to_locus[t_idx]] + 0.2
                link_labels = labels[q_idx]

                unique_groups = defaultdict(list)
                for e in range(len(q_idx)):
                    unique_groups[(link_labels[e], binned_idxs[e])].append(e)

                for (cluster_id, bin_idx), e_idxs in unique_groups.items():
                    color = self._get_cluster_color(cluster_id)
                    rgb = self._hex_to_rgb(color)
                    opacity = 0.4 + (bin_idx * 0.15)
                    rgba = f"rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {opacity:.2f})"

                    n_e = len(e_idxs)
                    x_poly = np.full((n_e, 6), np.nan, dtype=np.float64)
                    y_poly = np.full((n_e, 6), np.nan, dtype=np.float64)
                    e_arr = np.array(e_idxs, dtype=np.int32)

                    x_poly[:, 0] = q_starts[e_arr]
                    x_poly[:, 1] = q_ends[e_arr]
                    x_poly[:, 2] = t_ends[e_arr]
                    x_poly[:, 3] = t_starts[e_arr]
                    x_poly[:, 4] = q_starts[e_arr]

                    y_poly[:, 0] = q_y[e_arr]
                    y_poly[:, 1] = q_y[e_arr]
                    y_poly[:, 2] = t_y[e_arr]
                    y_poly[:, 3] = t_y[e_arr]
                    y_poly[:, 4] = q_y[e_arr]

                    name = f"Cluster {cluster_id}" if cluster_id >= 0 else "Unique"

                    fig.add_trace(
                        go.Scatter(
                            x=x_poly.flatten().tolist(),
                            y=y_poly.flatten().tolist(),
                            mode="lines",
                            fill="toself",
                            fillcolor=rgba,
                            line=dict(color="rgba(0,0,0,0)", width=0),
                            showlegend=False,
                            legendgroup=name,
                            hoverinfo="none",
                        )
                    )

                    # 4.5 Add Invisible Hover Markers for Links
                    cx = (q_starts[e_arr] + q_ends[e_arr] + t_starts[e_arr] + t_ends[e_arr]) / 4.0
                    cy = (q_y[e_arr] + t_y[e_arr]) / 2.0

                    link_texts = [
                        (
                            f"<b>Query:</b> {comparisons.locus_names[global_to_locus[q_idx[i]]]} "
                            f"({starts[q_idx[i]]}-{ends[q_idx[i]]})<br>"
                            f"<b>Target:</b> {comparisons.locus_names[global_to_locus[t_idx[i]]]} "
                            f"({starts[t_idx[i]]}-{ends[t_idx[i]]})<br>"
                            f"<b>Identity:</b> {pidents[i]:.1f}%"
                        )
                        for i in e_idxs
                    ]

                    fig.add_trace(
                        go.Scatter(
                            x=cx.tolist(),
                            y=cy.tolist(),
                            mode="markers",
                            marker=dict(size=1, color="rgba(0,0,0,0)"),
                            text=link_texts,
                            hoverinfo="text",
                            hoverlabel=dict(bgcolor=color),
                            legendgroup=name,
                            showlegend=False,
                        )
                    )

        # 5. Render Vectorised Gene Arrows & Hover Tooltips via GeneGlyphPlotter
        base_y = y_positions[global_to_locus]
        states_arr = (
            comparisons.gene_states if comparisons.gene_states is not None else np.zeros(total_genes, dtype=np.int8)
        )

        cluster_hover_texts = [
            f"<b>Locus:</b> {comparisons.locus_names[global_to_locus[i]]}<br>"
            f"<b>Gene:</b> {comparisons.gene_names[i]}<br>"
            f"<b>Product:</b> {comparisons.gene_descriptions[i]}<br>"
            f"<b>Coordinates:</b> {starts[i]}-{ends[i]}<br>"
            f"<b>State:</b> {GeneStyleManager.STATE_NAMES.get(int(states_arr[i]), 'NORMAL')}<br>"
            f"<b>Cluster:</b> {labels[i] if labels[i] >= 0 else 'Unique'}"
            for i in range(total_genes)
        ]

        GeneGlyphPlotter.render(
            fig=fig,
            starts=starts,
            ends=ends,
            strands=strands,
            y_positions=base_y,
            cluster_keys=labels,
            states=states_arr,
            hover_texts=cluster_hover_texts,
            dark_mode=dark_mode,
        )

        self._apply_standard_layout(
            fig=fig,
            yaxis=dict(
                tickmode="array",
                tickvals=y_positions.tolist(),
                ticktext=comparisons.locus_names,
                showgrid=False,
                zeroline=False,
            ),
            xaxis=dict(showgrid=False, zeroline=False),
            dark_mode=dark_mode,
        )

        return fig
