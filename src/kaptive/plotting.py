r"""Interactive Plotly visualization components for Kaptive serotyping and comparison maps.

Provides abstract base plotting utilities and concrete visualization implementations for:
- Serotyping locus gene diagrams ([`SerotypingResultPlotter`][kaptive.plotting.SerotypingResultPlotter])
- Multi-locus comparative synteny maps ([`LocusComparisonPlotter`][kaptive.plotting.LocusComparisonPlotter])

Classes:
    [`BasePlotter`][kaptive.plotting.BasePlotter]: Abstract base class providing color palette and geometry helpers.
    [`SerotypingResultPlotter`][kaptive.plotting.SerotypingResultPlotter]: Plotter for single-locus serotyping results.
    [`LocusComparisonPlotter`][kaptive.plotting.LocusComparisonPlotter]: Plotter for multi-locus comparison graphs.
"""

import zlib
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Sequence
from typing import Any

import numba
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from kaptive.compare import LocusComparisons
from kaptive.core.interval import Intervals
from kaptive.serotyping import GeneState, SerotypingResult


# Classes --------------------------------------------------------------------------------------------------------------
class BasePlotter(ABC):
    r"""Abstract base plotter providing common styling and geometry utilities.

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
            idx = zlib.crc32(cluster_id_or_name.encode()) % len(cls.PALETTE)
        return cls.PALETTE[idx]

    @classmethod
    def _hex_to_rgb(cls, hex_col: str) -> tuple[int, int, int]:
        r"""Convert a hex color string into an RGB integer tuple.

        Args:
            hex_col (str): Hexadecimal color string (e.g. '#FF0000').

        Returns:
            tuple[int, int, int]: RGB color tuple.
        """
        hex_col = hex_col.lstrip('#')
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

    @abstractmethod
    def __call__(self, *args: Any, **kwargs: Any) -> go.Figure:
        r"""Abstract call interface for rendering a Plotly figure.

        Returns:
            go.Figure: Interactive Plotly figure.
        """
        pass


class SerotypingResultPlotter(BasePlotter):
    r"""Plotter for visualising single-locus serotyping results.

    Generates interactive Plotly diagrams illustrating gene hit locations, locus contig
    backbones, gene state annotations, and detailed hover metadata.
    """

    STATE_NAMES = {s.value: s.name for s in GeneState}

    @classmethod
    def _get_state_styles(cls, dark_mode: bool = False) -> dict:
        r"""Get line and opacity style mappings for gene states.

        Args:
            dark_mode (bool): Whether dark mode palette styles should be used. Defaults to False.

        Returns:
            dict: Dictionary mapping [`GeneState`][kaptive.serotyping.GeneState] values
                to Plotly line and opacity style dicts.
        """
        line_color = "white" if dark_mode else "black"
        return {
            GeneState.NORMAL.value: dict(
                opacity=1.0, line=dict(width=1, dash="solid", color=line_color)
            ),
            GeneState.NOVEL.value: dict(
                opacity=0.5, line=dict(width=1, dash="dot", color=line_color)
            ),
            GeneState.PARTIAL.value: dict(
                opacity=0.7, line=dict(width=1, dash="dash", color="red")
            ),
            GeneState.TRUNCATED.value: dict(
                opacity=0.7, line=dict(width=1, dash="dash", color="orange")
            ),
        }

    @classmethod
    def __call__(cls, result: SerotypingResult, dark_mode: bool = False) -> go.Figure:
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

        # Generate polygons for genes
        lengths = plot_ends - plot_starts
        head_lens = np.minimum(lengths * 0.3, 400.0)

        is_rect = ~((states == GeneState.NORMAL.value) | (states == GeneState.NOVEL.value))

        y_positions = np.zeros(N, dtype=np.float64)

        x_coords, y_coords = cls._generate_gene_coordinates(
            plot_starts, plot_ends, plot_strands, y_positions, head_lens, is_rect
        )

        # Group traces by cluster and state
        cluster_names = genes.cluster_names
        traces_data = defaultdict(list)
        for i in range(N):
            traces_data[(cluster_names[i], states[i])].append(i)

        fig = go.Figure()

        # Add contig backbone
        piece_x = []
        piece_y = []
        tick_vals = []
        tick_text = []

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

            tick_vals.append((p_s + p_e) / 2)
            tick_text.append(label)

        fig.add_trace(
            go.Scatter(
                x=piece_x,
                y=piece_y,
                mode="lines",
                line=dict(color="white" if dark_mode else "gray", width=4),
                hoverinfo="none",
                showlegend=False,
            )
        )
        seen_clusters = set()

        for (c_name, st), idxs in traces_data.items():
            x_vals = x_coords[idxs].flatten().tolist()
            y_vals = y_coords[idxs].flatten().tolist()

            cluster_name = c_name
            color = cls._get_cluster_color(c_name)

            styles = cls._get_state_styles(dark_mode)
            style = styles.get(st, dict())

            show_leg = cluster_name not in seen_clusters
            seen_clusters.add(cluster_name)

            fig.add_trace(
                go.Scatter(
                    x=x_vals,
                    y=y_vals,
                    mode="lines",
                    fill="toself",
                    fillcolor=color,
                    line=style["line"],
                    opacity=style["opacity"],
                    name=cluster_name,
                    legendgroup=cluster_name,
                    showlegend=show_leg,
                    hoverinfo="none",
                )
            )

        # Add hover markers
        # Pre-resolve string lookups in list comprehensions which are faster than loops
        hover_texts = [
            f"<b>Gene:</b> {g}<br>"
            f"<b>Cluster:</b> {c}<br>"
            f"<b>Location:</b> {start}-{end} ({strand})<br>"
            f"<b>Coverage:</b> {cov:.1f}%<br>"
            f"<b>Identity:</b> {ident:.1f}%<br>"
            f"<b>State:</b> {cls.STATE_NAMES[s]}<br>"
            f"<b>Product:</b> {d}"
            for g, c, start, end, strand, cov, ident, s, d in zip(
                np.char.decode(genes.gene_ids),
                np.char.decode(genes.cluster_names),
                genes.t_starts,
                genes.t_ends,
                ["+" if st > 0 else "-" for st in genes.strands],
                genes.coverages,
                identities,
                states,
                np.char.decode(genes.product_descriptions),
            )
        ]

        hover_colors = [cls._get_cluster_color(c) for c in np.char.decode(genes.cluster_names)]

        fig.add_trace(
            go.Scatter(
                x=((plot_starts + plot_ends) / 2).tolist(),
                y=np.zeros(N).tolist(),
                mode="markers",
                marker=dict(size=1, color="rgba(0,0,0,0)"),
                text=hover_texts,
                hoverinfo="text",
                hoverlabel=dict(bgcolor=hover_colors),
                showlegend=False,
            )
        )

        fig.update_layout(
            title=f"Serotype Result: {result.best_locus_name} ({result.phenotype})",
            xaxis=dict(
                tickmode="array",
                tickvals=tick_vals,
                ticktext=tick_text,
                showgrid=False,
                zeroline=False,
            ),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False, range=[-1, 1]),
            template="plotly_dark" if dark_mode else "plotly_white",
            paper_bgcolor="rgba(0,0,0,0)" if dark_mode else "white",
            plot_bgcolor="rgba(0,0,0,0)" if dark_mode else "white",
            hovermode="closest",
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
            mask = (comparisons.edges.query_locus_indices == i - 1) & (
                comparisons.edges.target_locus_indices == i
            )
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
                min_coords.append(
                    np.min(comparisons.gene_intervals.starts[l_start : l_start + l_len])
                    + offsets[i]
                )

        if min_coords:
            offsets -= np.min(min_coords)

        return offsets.astype(np.int32)

    @classmethod
    def __call__(
        cls,
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
                offsets = cls._calculate_alignment_offsets(comparisons)
            else:
                offsets = np.zeros(n_loci, dtype=np.int32)
        else:
            offsets = np.asarray(offsets, dtype=np.int32)

        total_genes = len(comparisons.gene_intervals)
        if total_genes == 0:
            return go.Figure()

        edges = comparisons.edges

        # 1. Build Graph and cluster (connected components)
        if len(edges) > 0:
            labels = cls._connected_components(
                total_genes, edges.global_query_indices, edges.global_target_indices
            )

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

        starts = comparisons.gene_intervals.starts + offsets[global_to_locus]
        ends = comparisons.gene_intervals.ends + offsets[global_to_locus]
        strands = comparisons.gene_intervals.strands

        y_positions = -np.arange(n_loci, dtype=np.float64) * 2.0

        # 3. Plot Backbones (horizontal lines)
        for i in range(n_loci):
            l_start = comparisons.locus_offsets[i]
            l_len = comparisons.locus_lengths[i]
            if l_len > 0:
                min_x = np.min(starts[l_start : l_start + l_len])
                max_x = np.max(ends[l_start : l_start + l_len])
                fig.add_trace(
                    go.Scatter(
                        x=[min_x, max_x],
                        y=[y_positions[i], y_positions[i]],
                        mode="lines",
                        line=dict(color="gray", width=2),
                        showlegend=False,
                        hoverinfo="none",
                    )
                )

        # 4. Add Homology Links (Polygons with opacity binning)
        if len(edges) > 0:
            mask = np.ones(len(edges), dtype=bool)
            if not show_all_links:
                mask = edges.target_locus_indices == edges.query_locus_indices + 1

            q_idx = edges.global_query_indices[mask]
            t_idx = edges.global_target_indices[mask]
            pidents = edges.alignments.pidents[mask]

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
                    color = cls._get_cluster_color(cluster_id)
                    rgb = cls._hex_to_rgb(color)
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

        # 5. Add Vectorised Gene Arrows
        lengths = ends - starts
        head_lens = np.minimum(lengths * 0.3, 400.0)

        base_y = y_positions[global_to_locus]

        x_coords, y_coords = cls._generate_gene_coordinates(
            starts, ends, strands, base_y, head_lens, is_rect=None
        )

        # Group arrows by cluster trace
        traces_data = defaultdict(list)
        for i in range(total_genes):
            traces_data[labels[i]].append(i)

        for cluster_id, idxs in traces_data.items():
            x_vals = x_coords[idxs].flatten().tolist()
            y_vals = y_coords[idxs].flatten().tolist()

            color = cls._get_cluster_color(cluster_id)
            name = f"Cluster {cluster_id}" if cluster_id >= 0 else "Unique"

            fig.add_trace(
                go.Scatter(
                    x=x_vals,
                    y=y_vals,
                    mode="lines",
                    fill="toself",
                    fillcolor=color,
                    line=dict(color="white" if dark_mode else "black", width=1),
                    name=name,
                    legendgroup=name,
                    showlegend=True,
                    hoverinfo="none",
                )
            )

            # 6. Invisible Hover Markers
            cx = (starts[idxs] + ends[idxs]) / 2.0
            cy = base_y[idxs]

            cluster_hover_texts = [
                f"<b>Locus:</b> {comparisons.locus_names[global_to_locus[i]]}<br>"
                f"<b>Gene:</b> {comparisons.gene_names[i]}<br>"
                f"<b>Description:</b> {comparisons.gene_descriptions[i]}<br>"
                f"<b>Coordinates:</b> {starts[i]}-{ends[i]}<br>"
                f"<b>Cluster:</b> {cluster_id if cluster_id >= 0 else 'Unique'}"
                for i in idxs
            ]

            fig.add_trace(
                go.Scatter(
                    x=cx.tolist(),
                    y=cy.tolist(),
                    mode="markers",
                    marker=dict(size=1, color="rgba(0,0,0,0)"),
                    text=cluster_hover_texts,
                    hoverinfo="text",
                    hoverlabel=dict(bgcolor=color),
                    legendgroup=name,
                    showlegend=False,
                )
            )

        fig.update_layout(
            yaxis=dict(
                tickmode="array",
                tickvals=y_positions.tolist(),
                ticktext=comparisons.locus_names,
                showgrid=False,
                zeroline=False,
            ),
            xaxis=dict(showgrid=False, zeroline=False),
            template="plotly_dark" if dark_mode else "plotly_white",
            paper_bgcolor="rgba(0,0,0,0)" if dark_mode else "white",
            plot_bgcolor="rgba(0,0,0,0)" if dark_mode else "white",
            hovermode="closest",
        )

        return fig

