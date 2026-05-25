from typing import Iterable, NamedTuple
from collections import defaultdict

import numpy as np

from kaptive.core.interval import Strand
from kaptive.core.alignment import AlignmentBatch, AlignmentRecord


# Classes --------------------------------------------------------------------------------------------------------------
class Edge(NamedTuple):
    """
    Represents a directed link between two contigs in an assembly graph.

    Attributes:
        u: Source contig name.
        u_strand: Strand of the source contig (Strand.FORWARD or Strand.REVERSE).
        v: Destination contig name.
        v_strand: Strand of the destination contig.
        overlap: Length of the sequence overlap between the two contigs in base pairs.
    """
    u: str
    u_strand: Strand
    v: str
    v_strand: Strand
    overlap: int = 0

    def reverse(self) -> 'Edge':
        """Returns a new Edge object representing the reverse traversal (v -> u)."""
        return Edge(self.v, self.v_strand, self.u, self.u_strand, self.overlap)


class Graph:
    """
    Manages a graph, supporting both directed and undirected traversal.

    Nodes represent contigs (or fragments), and edges represent physical
    adjacencies (e.g., from an assembly GFA file).
    """
    __slots__ = ('adj', 'in_adj', 'edges', '_nodes', 'directed')
    def __init__(self, edges: Iterable[Edge] | None = None, directed: bool = True):
        """
        Initialize the graph.

        Args:
            edges: Optional iterable of Edge objects to seed the graph.
            directed: If True, edges are strictly one-way. If False, every
                added edge u -> v implicitly adds v -> u.
        """
        # Adjacency list: maps node ID to a set of outgoing Edge objects *starting* uom that node.
        # For undirected graphs, this will include edges representing reverse traversal.
        self.adj: dict[str, set[Edge]] = defaultdict(set)
        # In-degree adjacency list for efficient reverse lookups
        self.in_adj: dict[str, set[Edge]] = defaultdict(set)
        # Set of unique Edge objects fundamentally added to the graph.
        self.edges: set[Edge] = set()
        self._nodes: set[str] = set()
        self.directed: bool = directed
        if edges is not None:
            for edge in edges:
                self.add_edge(edge)

    def __repr__(self):
        # Note: len(self.edges) counts only the *unique* edge objects added,
        # not the total number of traversable connections in the undirected case.
        return f"{'Directed' if self.directed else 'Undirected'} Graph with {len(self._nodes)} nodes and {len(self.edges)} defined edges"

    def __iter__(self):
        return iter(self.edges)

    def __len__(self):
        return len(self.edges)

    def add_node(self, node: str):
        """Adds a node to the graph if it doesn't already exist."""
        self._nodes.add(node)

    def add_edge(self, edge: Edge):
        """
        Adds an edge to the graph.

        If the graph is undirected, the reverse connection is also added
        to the adjacency lists to allow bidirectional traversal.
        """
        # Add the nodes to the set of known node IDs
        self.add_node(edge.u)
        self.add_node(edge.v)

        # This helps track the originally added edges vs implicit reverse ones
        self.edges.add(edge)  # Add the primary edge representation if it's new

        # Check if this specific edge object is already in the adjacency list for the 'u' node
        self.adj[edge.u].add(edge)
        self.in_adj[edge.v].add(edge)

        # If the graph is undirected, add reverse connectivity as well
        if not self.directed:
            # Create a conceptual reverse edge for traversal and add to the adjacency list of the 'v' node
            reverse_edge = edge.reverse()
            self.in_adj[edge.u].add(reverse_edge)
            self.adj[edge.v].add(reverse_edge)
            # Note: We do not add the reverse_edge to self.edges unless it's explicitly added later by the user

    def get_neighbors(self, node_id: str) -> set[Edge]:
        """
        Returns the set of outgoing edges for a given node ID.

        Respects graph directionality. For undirected graphs, this includes
        implicit reverse connections.
        """
        return self.adj.get(node_id, set())


class AssemblyGraph(Graph):
    __slots__ = ('_contigs',)
    def __init__(self, genome: 'GenomeAssembly'):
        super().__init__(genome.edges, directed=False)
        self._contigs = genome.contigs

    @property
    def contigs(self):
        return self._contigs

    def find_bounded_paths(self, start_ctg: str, start_strand: Strand, target_ctg: str,
                           target_strand: Strand, expected_len: int, tolerance: int) -> list[dict]:
        """
        Finds all physical paths between two contigs within a length constraint.

        Args:
            start_ctg: Starting contig ID.
            start_strand: Strand to exit the start contig from.
            target_ctg: Target contig ID.
            target_strand: Required strand to enter the target contig.
            expected_len: The gap distance observed in the query sequence.
            tolerance: bp tolerance for the path length matching the expected length.

        Returns:
            list[dict]: Valid paths found, with 'contigs', 'length', and 'min_depth'.
        """
        # Stack payload: (current_contig, exit_strand, path_list, accumulated_len, bottleneck_depth)
        stack = [(start_ctg, start_strand, [start_ctg], 0, float('inf'))]
        valid_paths = []

        while stack:
            curr_ctg, curr_strand, path, dist, min_dp = stack.pop()

            # Base Case: Reached the Sink anchor
            if curr_ctg == target_ctg:
                if curr_strand == target_strand:  # Did we arrive on the correct biological strand?
                    valid_paths.append({'contigs': path, 'length': dist, 'min_depth': min_dp})
                continue

            # Prune: We have wandered too far down a dead end
            if dist > expected_len + tolerance:
                continue

            # Graph Traversal
            for edge in self.get_neighbors(curr_ctg):
                if edge.u_strand != curr_strand:
                    continue

                n_ctg = edge.v
                if n_ctg in path:
                    continue  # Prevent cyclic infinite loops

                record = self._contigs[n_ctg]
                n_len = record.length if record else 0
                n_dp = record.depth if record else 1.0
                overlap_len = getattr(edge, 'overlap', 0)

                # CRITICAL FIX: If this is the target anchor, its length doesn't belong in the gap!
                added_dist = (n_len - overlap_len) if n_ctg != target_ctg else -overlap_len

                stack.append((
                    n_ctg,
                    edge.v_strand,
                    path + [n_ctg],
                    dist + added_dist,
                    min(min_dp, n_dp)
                ))

        return valid_paths

    def stitch_alignments(self, alignments: 'AlignmentBatch') -> 'AlignmentBatch':
        """
        Stitch together partial alignments that span multiple graph nodes.
        Takes a single concatenated SoA batch and returns an updated SoA batch.
        """
        # 1. Vectorized identification of partials (Instant)
        partial_mask = alignments.is_partial_mask
        partial_indices = np.where(partial_mask)[0]

        if len(partial_indices) == 0:
            return alignments  # Nothing to stitch!

        # 2. Extract ONLY the partials into Records for the DFS
        partial_alns = defaultdict(list)
        for idx in partial_indices:
            rec = alignments.get_record(idx)
            partial_alns[rec.q_name].append(rec)

        resolved_paths = []
        used_indices = set()  # Track the original SoA indices consumed by the stitcher

        # 3. Fragment Chaining via DAG DFS (Logic remains mostly the same)
        for q_name, fragments in partial_alns.items():
            fragments.sort(key=lambda x: x.q_start)
            used_in_this_qname = set()

            for i in range(len(fragments)):
                if i in used_in_this_qname: continue
                start_frag = fragments[i]
                stack = [(i, [start_frag], {i})]

                best_chain = []
                best_chain_len = 0
                best_chain_used = set()

                while stack:
                    curr_idx, curr_path, curr_used = stack.pop()
                    curr_frag = curr_path[-1]
                    extended = False

                    for j in range(len(fragments)):
                        if j in curr_used or j in used_in_this_qname: continue
                        next_frag = fragments[j]
                        expected_gap = next_frag.q_start - curr_frag.q_end

                        if expected_gap < -50: continue

                        paths = self.find_bounded_paths(
                            curr_frag.t_name, curr_frag.strand,
                            next_frag.t_name, next_frag.strand,
                            expected_gap, tolerance=2000
                        )

                        if paths:
                            # ... (Your existing scoring logic to find best_p_idx) ...
                            best_p_idx = 0  # Placeholder for brevity
                            winning_contigs = paths[best_p_idx]['contigs']
                            extension = self._build_stitching_payload(curr_frag, next_frag, winning_contigs)

                            new_path = curr_path + extension[1:]
                            stack.append((j, new_path, curr_used | {j}))
                            extended = True

                    if not extended:
                        chain_cov = curr_path[-1].q_end - curr_path[0].q_start
                        if chain_cov > best_chain_len:
                            best_chain = curr_path
                            best_chain_len = chain_cov
                            best_chain_used = curr_used

                if len(best_chain_used) > 1:
                    resolved_paths.append(best_chain)
                    used_in_this_qname.update(best_chain_used)
                    for f in best_chain:
                        if f.idx != -1:
                            used_indices.add(f.idx)

        # 4. Rebuild the Master Batch
        if not resolved_paths:
            return alignments

        # Create a boolean mask keeping everything EXCEPT the stitched fragments
        keep_mask = np.ones(len(alignments), dtype=bool)
        if used_indices:
            keep_mask[list(used_indices)] = False

        intact_batch = alignments.filter(keep_mask)

        # Convert the newly stitched paths back into an AlignmentBatch
        # We flatten the list of lists into a single list of synthetic AlignmentRecords
        stitched_records = [rec for path in resolved_paths for rec in path]

        if stitched_records:
            # You'll need a helper classmethod to convert a list of Records back to a Batch
            stitched_batch = AlignmentBatch.from_records(stitched_records)
            return AlignmentBatch.concat([intact_batch, stitched_batch])

        return intact_batch

    def _build_stitching_payload(self, h_u: 'AlignmentRecord', h_v: 'AlignmentRecord', path_contigs: list[str]) -> list[
        'AlignmentRecord']:
        """
        Converts a list of graph contig names into a continuous sequence of AlignmentRecords.

        Generates 'synthetic' records for unaligned intermediate nodes so they
        can be processed by the LocusBuilder as part of a single path.
        """
        payload = [h_u]

        # Iterate over only the intermediate contigs (excluding the h_u and h_v anchors)
        for ctg in path_contigs[1:-1]:
            record = self.contigs.get(ctg)
            ctg_len = record.length if record else 0

            # Create a synthetic alignment that claims the entire unaligned contig
            payload.append(AlignmentRecord(
                idx=-1,  # Flag as a synthetic/mock record
                q_name=h_u.q_name,
                q_length=h_u.q_length,
                q_start=h_u.q_end,  # Conceptually sits between the anchors
                q_end=h_v.q_start,
                t_name=ctg,
                t_length=ctg_len,
                t_start=0,  # The entire contig is part of the path
                t_end=ctg_len,
                strand=Strand.FORWARD,  # Default to forward for the traversal sequence
                length=ctg_len,
                match=0,
                mismatch=0,
                quality=0,
                cigar="*"  # No CIGAR exists for synthetic nodes
            ))

        payload.append(h_v)

        return payload
