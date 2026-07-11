from kaptive.core.interval import Intervals
from kaptive.compare import LocusComparisonEdges, LocusComparisons
from kaptive.core.pairwise import PairwiseAlignments
from kaptive.plotting import LocusComparisonPlotter
import numpy as np

# dummy
b1 = Intervals.from_dict({'starts': [0, 1000], 'ends': [500, 1500], 'strands': [1, -1]})
b2 = Intervals.from_dict({'starts': [100, 1100], 'ends': [600, 1600], 'strands': [1, -1]})

alignments = PairwiseAlignments(
    scores=np.array([100, 100]),
    matches=np.array([500, 500]),
    mismatches=np.array([0, 0]),
    gaps=np.array([0, 0]),
    q_starts=np.array([0, 0]),
    q_ends=np.array([500, 500]),
    t_starts=np.array([0, 0]),
    t_ends=np.array([500, 500])
)

edges = LocusComparisonEdges(
    query_locus_indices=np.array([0, 0]),
    target_locus_indices=np.array([1, 1]),
    query_indices=np.array([0, 1]),
    target_indices=np.array([0, 1]),
    global_query_indices=np.array([0, 1]),
    global_target_indices=np.array([2, 3]),
    alignments=alignments
)

comparisons = LocusComparisons(
    edges=edges,
    locus_names=("Locus 1", "Locus 2"),
    locus_lengths=np.array([2, 2], dtype=np.int32),
    locus_offsets=np.array([0, 2], dtype=np.int32),
    gene_names=np.array(["A", "B", "C", "D"], dtype=object),
    gene_descriptions=np.array(["", "", "", ""], dtype=object),
    gene_intervals=Intervals(
        starts=np.concatenate([b1.starts, b2.starts]),
        ends=np.concatenate([b1.ends, b2.ends]),
        strands=np.concatenate([b1.strands, b2.strands]),
        original_indices=np.arange(4, dtype=np.int32)
    )
)

plotter = LocusComparisonPlotter()
fig = plotter(comparisons)
print("Success")
