"""
Module to handle query (contigs) and target (features) IO.
"""
from dataclasses import dataclass, field
from collections.abc import Iterator
from pathlib import Path
from re import compile as re_compile
from typing import IO, ClassVar, Callable,  Iterable, NamedTuple
from collections import defaultdict

# NOTE: When we enforce python>=3.14 we can use the compression library!
from bz2 import open as bzopen
from gzip import open as gzopen
from lzma import open as lzopen

import numpy as np
import numpy.typing as npt

from kaptive.core.seq import SeqRecord, Sequences
from kaptive.core.interval import Strand
from kaptive.core.alignment import Alignments, Alignment


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class GenomeAssembly:
    """
    Container for a genome assembly, including contigs and their graph topology.

    Handles FASTA and GFA formats, with support for transparent decompression.
    """
    _SEQUENCE_FILE_REGEX = re_compile(
        r'\.('
        r'(?P<fasta>f(asta|a|na|fn|as|aa))|'
        r'(?P<gfa>gfa)|'
        r')\.?(?P<compression>(gz|bz2|xz))?$'
    )
    _OPENERS: ClassVar[dict[str, Callable]] = {'gz': gzopen, 'bz2': bzopen, 'xz': lzopen}
    id: str
    contigs: Sequences
    edges: Edges
    id_map: dict[str, int] = field(init=False, repr=False, hash=False, compare=False)
    node_bounds: npt.NDArray[np.int32] = field(init=False, repr=False, hash=False, compare=False)

    def __post_init__(self):
        object.__setattr__(self, 'id_map', {name: i for i, name in enumerate(self.contigs.ids)})

        if len(self.edges) == 0:
            bounds = np.zeros((len(self.contigs), 2), dtype=np.int32)
            bounds.flags.writeable = False
            object.__setattr__(self, 'node_bounds', bounds)
        else:
            # Assembly graphs are physically undirected; create bidirectional edges
            all_edges = Edges.concat([self.edges, self.edges.reverse()])

            # Map source names to integer indices for rapid sorting
            u_ints = np.array([self.id_map[u] for u in all_edges.u_names], dtype=np.int32)

            # Lexsort to group by source node
            order = np.argsort(u_ints)
            sorted_edges = all_edges[order]
            sorted_u = u_ints[order]

            # Precompute adjacency slices for O(1) neighbor lookups
            bounds = np.zeros((len(self.contigs), 2), dtype=np.int32)
            unique_u, indices = np.unique(sorted_u, return_index=True)
            bounds[unique_u, 0] = indices
            bounds[unique_u, 1] = np.append(indices[1:], len(sorted_edges))

            bounds.flags.writeable = False
            object.__setattr__(self, 'edges', sorted_edges)
            object.__setattr__(self, 'node_bounds', bounds)

    @classmethod
    def ensure(cls, genome: GenomeAssembly | str | Path) -> GenomeAssembly:
        """Ensures the input is a GenomeAssembly, loading it from file if necessary.
        
        Args:
            genome: A GenomeAssembly instance, or a path to a FASTA/GFA file.
            
        Returns:
            GenomeAssembly: The validated or loaded assembly.
        """
        if isinstance(genome, cls):
            return genome
        return cls.from_file(genome)

    def __len__(self) -> int:
        """Total number of base pairs in the assembly."""
        return len(self.contigs.seqs)

    def __iter__(self) -> Iterator[SeqRecord]:
        return iter(self.contigs)

    def __str__(self):
        return self.id

    def __getitem__(self, item: str) -> bytes:
        """Access a contig sequence by its ID."""
        idx = self.id_map[item]
        s = self.contigs.offsets[idx]
        l = self.contigs.lengths[idx]
        return self.contigs.seqs[s:s + l].tobytes()

    @classmethod
    def from_file(cls, file: str | Path):
        """
        Load an assembly from a FASTA or GFA file.

        Args:
            file: Path to the file. Supports .gz, .bz2, and .xz compression.
        """
        file = Path(file) # type: Path
        if not (m := cls._SEQUENCE_FILE_REGEX.search(file.name)):
            raise NotImplementedError(f'Unsupported format: {file}')
        
        if m.group('fasta'):
            from kaptive.io import FastaReader as reader
        else:
            from kaptive.io import GfaReader as reader
            
        with cls._OPENERS.get(m.group('compression'), open)(file, mode='rb') as handle:
            return cls.from_stream(handle, reader, file.name.rstrip(m.group()))

    @classmethod
    def from_stream(cls, handle: IO[bytes], reader, id_: str | None = None):
        """Load an assembly from an open file stream using the specified reader."""
        contigs, edges = [], []
        for record in reader(handle):
            if isinstance(record, SeqRecord):
                contigs.append(record)
            elif isinstance(record, Edge):
                edges.append(record)
        return cls(id_ or handle.name, Sequences.from_records(contigs), Edges.from_edges(edges))

    def get_neighbors(self, node: str | int) -> Edges:
        """Returns all outgoing edges for a given contig in O(1) time as an Edges slice."""
        node_idx = self.id_map[node] if isinstance(node, str) else node
        s, e = self.node_bounds[node_idx]
        if s == e:
            return Edges.empty()
        return self.edges[s:e]

    def find_bounded_paths(self, start_ctg: str, start_strand: Strand, target_ctg: str,
                           target_strand: Strand, expected_len: int, tolerance: int) -> Paths:
        """Finds all simple paths between two contigs that satisfy a length constraint.

        This method performs a depth-first search (DFS) to explore all possible paths of contigs
        connecting a `start_ctg` to a `target_ctg`. The search is constrained by the `expected_len`
        (the gap distance observed between two alignments on a query sequence) and a `tolerance`.
        It prunes paths that become too long and ensures the final path arrives at the target contig
        on the correct strand.

        Args:
            start_ctg (str): The identifier of the starting contig.
            start_strand (Strand): The strand to exit from the starting contig.
            target_ctg (str): The identifier of the target contig.
            target_strand (Strand): The strand required to enter the target contig.
            expected_len (int): The expected length of the path (sum of intermediate contig lengths
                minus overlaps).
            tolerance (int): The allowable deviation from `expected_len`.

        Returns:
            Paths: A high-performance batch of all valid paths found.
        """
        # Convert names to integer indices to avoid slow string comparisons during traversal
        start_idx = self.id_map[start_ctg]
        target_idx = self.id_map[target_ctg]

        # Stack payload: (current_contig_idx, exit_strand, path_list, accumulated_len, bottleneck_depth)
        stack = [(start_idx, start_strand, [start_idx], 0, float('inf'))]
        valid_nodes, valid_dists, valid_depths = [], [], []

        lengths = self.contigs.lengths
        has_depths = hasattr(self.contigs, 'depths')

        while stack:
            curr_ctg, curr_strand, path, dist, min_dp = stack.pop()

            # Base Case: Reached the Sink anchor
            if curr_ctg == target_idx:
                if curr_strand == target_strand:  # Did we arrive on the correct biological strand?
                    valid_nodes.append(path)
                    valid_dists.append(dist)
                    valid_depths.append(min_dp)
                continue

            # Prune: We have wandered too far down a dead end
            if dist > expected_len + tolerance:
                continue

            # Fetch O(1) neighbors
            neighbors = self.get_neighbors(curr_ctg)
            if len(neighbors) == 0:
                continue

            # Vectorized filtering for the correct outgoing strand
            valid_mask = neighbors.u_strands == curr_strand.value
            if not np.any(valid_mask):
                continue

            valid_neighbors = neighbors[valid_mask]

            # Graph Traversal
            for i in range(len(valid_neighbors)):
                n_ctg_str = valid_neighbors.v_names[i]
                n_ctg = self.id_map[n_ctg_str]

                if n_ctg in path:
                    continue  # Prevent cyclic infinite loops

                v_strand = Strand(valid_neighbors.v_strands[i])
                overlap_len = valid_neighbors.overlaps[i]

                n_len = lengths[n_ctg]
                n_dp = self.contigs.depths[n_ctg] if has_depths else 1.0

                # CRITICAL FIX: If this is the target anchor, its length doesn't belong in the gap!
                added_dist = (n_len - overlap_len) if n_ctg != target_idx else -overlap_len

                stack.append((
                    n_ctg,
                    path + [n_ctg],
                    dist + added_dist,
                    min(min_dp, n_dp)
                ))

        return Paths.from_lists(valid_nodes, valid_dists, valid_depths)

    def stitch_alignments(self, alignments: Alignments) -> Alignments:
        """Stitches together fragmented alignments that span multiple contigs in the graph.

        This method identifies alignments that are partial (i.e., do not cover the full query sequence)
        and attempts to connect them by finding valid paths in the assembly graph between the contigs
        they serotyping to. If a valid path is found, it creates a new, "stitched" alignment record
        representing the complete path.

        The process is as follows:
        
        1.  Identifies all partial alignments in the input `Alignments`.
        2.  Groups these partials by their query name.
        3.  For each group, it attempts to chain fragments together by using `find_bounded_paths`
            to bridge the gaps between them.
        4.  If a successful chain is found, the original partial alignments are removed, and a new
            set of records (including synthetic ones for intermediate contigs) is created.
        5.  Returns a new `Alignments` containing the original complete alignments plus the
            newly stitched alignments.

        Args:
            alignments (Alignments): A batch of alignments to be processed.

        Returns:
            Alignments: An updated batch with fragmented alignments stitched into longer ones.
        """
        # 1. Vectorized identification of partials (Instant)
        partial_mask = alignments.is_partial_mask
        partial_indices = np.where(partial_mask)[0]

        if len(partial_indices) == 0:
            return alignments  # Nothing to stitch!

        # 2. Extract ONLY the partials into Records for the DFS
        partial_alns = defaultdict(list)
        for idx in partial_indices:
            rec = alignments[idx]
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

                        if len(paths) > 0:
                            # Select the path whose cumulative length best matches the expected physical gap size
                            best_p_idx = int(np.argmin(np.abs(paths.distances - expected_gap)))
                            
                            s, c = paths.offsets[best_p_idx], paths.counts[best_p_idx]
                            winning_nodes = paths.nodes[s:s + c]
                            extension = self._build_stitching_payload(curr_frag, next_frag, winning_nodes)

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

        # Convert the newly stitched paths back into an Alignments
        # We flatten the list of lists into a single list of synthetic AlignmentRecords
        stitched_records = [rec for path in resolved_paths for rec in path]

        if stitched_records:
            # You'll need a helper classmethod to convert a list of Records back to a Batch
            stitched_batch = Alignments.from_records(stitched_records)
            return Alignments.concat([intact_batch, stitched_batch])

        return intact_batch

    def _build_stitching_payload(self, h_u: Alignment, h_v: Alignment, path_nodes: npt.NDArray[np.int32]
                                 ) -> list[Alignment]:
        """Converts a list of contig names from a graph path into a continuous list of `AlignmentRecord`s.

        This helper method takes the start and end alignment fragments (`h_u`, `h_v`) and the list of
        contigs that form a path between them. It generates "synthetic" `AlignmentRecord`s for the
        intermediate, unaligned contigs in the path. These synthetic records serve as placeholders,
        allowing the entire path to be treated as a single, continuous alignment by downstream logic.

        Args:
            h_u (Alignment): The starting alignment fragment.
            h_v (Alignment): The ending alignment fragment.
            path_nodes (npt.NDArray[np.int32]): The array of contig indices forming the path.

        Returns:
            list[AlignmentRecord]: A list of records representing the full stitched path, including
                the original start and end fragments and synthetic records for the gaps.
        """
        payload = [h_u]

        # Iterate over only the intermediate contigs (excluding the h_u and h_v anchors)
        for ctg_idx in path_nodes[1:-1]:
            ctg_name = self.contigs.ids[ctg_idx]
            ctg_len = self.contigs.lengths[ctg_idx]

            # Create a synthetic alignment that claims the entire unaligned contig
            payload.append(Alignment(
                idx=-1,  # Flag as a synthetic/mock record
                q_name=h_u.q_name,
                q_length=h_u.q_length,
                q_start=h_u.q_end,  # Conceptually sits between the anchors
                q_end=h_v.q_start,
                t_name=ctg_name,
                t_length=ctg_len,
                t_start=0,  # The entire contig is part of the path
                t_end=ctg_len,
                strand=Strand.FORWARD,  # Default to forward for the traversal sequence
                length=ctg_len,
                match=0,
                mismatch=0,
                quality=0,
                cigar=np.empty(0, dtype=np.uint32)  # No CIGAR exists for synthetic nodes
            ))

        payload.append(h_v)

        return payload


@dataclass(frozen=True, slots=True)
class Paths:
    """A high-performance SoA container for assembly graph paths.

    Stores a collection of paths (each being a sequence of contig indices) in a flat memory layout,
    along with their corresponding lengths and bottleneck depths.

    Attributes:
        nodes (npt.NDArray[np.int32]): A flat array of all contig indices across all paths.
        offsets (npt.NDArray[np.int32]): The starting index in `nodes` for each path.
        counts (npt.NDArray[np.int32]): The number of contigs in each path.
        distances (npt.NDArray[np.int32]): The cumulative physical length of each path.
        min_depths (npt.NDArray[np.float64]): The minimum (bottleneck) read depth of each path.
    """
    nodes: npt.NDArray[np.int32]
    offsets: npt.NDArray[np.int32]
    counts: npt.NDArray[np.int32]
    distances: npt.NDArray[np.int32]
    min_depths: npt.NDArray[np.float64]

    def __len__(self) -> int:
        return len(self.distances)

    @classmethod
    def empty(cls) -> 'Paths':
        return cls(
            np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.float64)
        )

    @classmethod
    def from_lists(cls, nodes: list[list[int]], distances: list[int], depths: list[float]) -> 'Paths':
        if not nodes:
            return cls.empty()
        
        counts = np.array([len(p) for p in nodes], dtype=np.int32)
        offsets = np.zeros(len(counts), dtype=np.int32)
        if len(counts) > 1:
            np.cumsum(counts[:-1], out=offsets[1:])
            
        flat_nodes = np.concatenate([np.array(p, dtype=np.int32) for p in nodes])
        
        return cls(
            nodes=flat_nodes, offsets=offsets, counts=counts,
            distances=np.array(distances, dtype=np.int32),
            min_depths=np.array(depths, dtype=np.float64)
        )


class Edge(NamedTuple):
    """A lightweight, immutable representation of a directed link between two nodes in a graph.

    In the context of assembly graphs, nodes are contigs. This `NamedTuple` captures the essential
    information about a physical adjacency between two contigs, including the orientation (strand)
    at which the connection occurs and the length of the sequence overlap.

    Attributes:
        u (str): The name of the source contig (node).
        u_strand (Strand): The strand of the source contig from which the edge originates
            (e.g., `Strand.FORWARD` means the edge comes off the 3' end of the forward strand).
        v (str): The name of the destination contig (node).
        v_strand (Strand): The strand of the destination contig at which the edge arrives.
        overlap (int): The length of the sequence overlap between the two contigs in base pairs.
    """
    u: str
    u_strand: Strand
    v: str
    v_strand: Strand
    overlap: int = 0

    def reverse(self) -> Edge:
        """Returns a new Edge object representing the reverse traversal of this edge (v -> u).

        This is useful for navigating undirected graphs where every connection is bidirectional.

        Returns:
            Edge: A new `Edge` with the source and destination nodes and strands swapped.
        """
        return Edge(self.v, self.v_strand, self.u, self.u_strand, self.overlap)


@dataclass(frozen=True, slots=True)
class Edges:
    """A high-performance SoA container for assembly graph edges.

    This class stores a collection of edges in a Structure-of-Arrays (SoA) layout using NumPy arrays.
    It provides significant performance benefits over a list of `Edge` objects, especially when
    filtering or performing vectorized operations on large graphs.

    Attributes:
        u_names (npt.NDArray[np.object_]): Array of source contig names.
        u_strands (npt.NDArray[np.int8]): Array of source contig strands (1 for FORWARD, -1 for REVERSE).
        v_names (npt.NDArray[np.object_]): Array of destination contig names.
        v_strands (npt.NDArray[np.int8]): Array of destination contig strands.
        overlaps (npt.NDArray[np.int32]): Array of overlap lengths in base pairs.
    """
    u_names: npt.NDArray[np.object_]
    u_strands: npt.NDArray[np.int8]
    v_names: npt.NDArray[np.object_]
    v_strands: npt.NDArray[np.int8]
    overlaps: npt.NDArray[np.int32]

    def __len__(self) -> int:
        """Returns the number of edges in the batch."""
        return len(self.u_names)

    @classmethod
    def empty(cls) -> Edges:
        """Creates an empty Edges collection."""
        return cls(
            np.empty(0, dtype=object),
            np.empty(0, dtype=np.int8),
            np.empty(0, dtype=object),
            np.empty(0, dtype=np.int8),
            np.empty(0, dtype=np.int32)
        )

    @classmethod
    def from_edges(cls, edges: Iterable[Edge]) -> Edges:
        """Constructs an Edges object from an iterable of Edge objects."""
        # OPTIMIZATION: Fast C-level list comprehension + zip extraction
        data = [(e.u, e.u_strand, e.v, e.v_strand, e.overlap) for e in edges]
        if not data:
            return cls.empty()

        u, us, v, vs, o = zip(*data, strict=True)
        return cls(
            np.array(u, dtype=object),
            np.array(us, dtype=np.int8),
            np.array(v, dtype=object),
            np.array(vs, dtype=np.int8),
            np.array(o, dtype=np.int32)
        )

    @classmethod
    def concat(cls, batches: Iterable[Edges]) -> Edges:
        """Concatenates multiple Edges objects into a single, larger collection."""
        batches = list(batches)
        if not batches:
            raise ValueError("Cannot concatenate an empty iterable of batches")

        return cls(
            np.concatenate([b.u_names for b in batches]),
            np.concatenate([b.u_strands for b in batches]),
            np.concatenate([b.v_names for b in batches]),
            np.concatenate([b.v_strands for b in batches]),
            np.concatenate([b.overlaps for b in batches])
        )

    def __getitem__(self, item) -> Edge | Edges:
        """Accesses edges by index, slice, or boolean mask.

        - If `item` is an integer, it returns a single `Edge` object.
        - If `item` is a slice or boolean mask, it returns a new, smaller `Edges` collection.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return Edge(
                u=self.u_names[item],
                u_strand=Strand(self.u_strands[item]),
                v=self.v_names[item],
                v_strand=Strand(self.v_strands[item]),
                overlap=self.overlaps[item]
            )

        return Edges(
            self.u_names[item],
            self.u_strands[item],
            self.v_names[item],
            self.v_strands[item],
            self.overlaps[item]
        )

    def filter(self, mask: np.ndarray) -> Edges:
        """Returns a new batch containing only elements where mask is True. (Alias for self[mask])"""
        return self[mask]

    def reverse(self) -> Edges:
        """Returns a new Edges collection representing the reverse traversal of all edges (v -> u)."""
        return Edges(
            u_names=self.v_names,
            u_strands=self.v_strands,
            v_names=self.u_names,
            v_strands=self.u_strands,
            overlaps=self.overlaps
        )