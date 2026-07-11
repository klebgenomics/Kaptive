import json
import logging
from kaptive.serotyping import SerotypingResult, GeneState, GeneHits, LocusPieces
from kaptive.core.seq import Sequences
from kaptive.plotting import SerotypingResultPlotter
import numpy as np

def main():
    hits = GeneHits(
        t_indices=np.array([0, 0]),
        q_indices=np.array([0, 1]),
        t_starts=np.array([100, 500]),
        t_ends=np.array([300, 700]),
        q_starts=np.array([0, 0]),
        q_ends=np.array([200, 200]),
        strands=np.array([1, -1]),
        q_lens=np.array([200, 200]),
        gene_ids=["gene1", "gene2"],
        cluster_names=["A", "B"],
        expected_positions=np.array([1.0, 2.0]),
        expected_strands=np.array([1, 1]),
        is_expected=np.array([True, True]),
        is_extra=np.array([False, False]),
        is_inside=np.array([True, True]),
    )
    
    pieces = LocusPieces(
        ctg_indices=np.array([0]),
        starts=np.array([0]),
        ends=np.array([1000]),
        strands=np.array([1])
    )
    
    seqs = Sequences(ids=["ctg1"], seqs=["A"*1000], descriptions=[""])
    
    res = SerotypingResult(
        genome="genome1",
        best_locus_name="locus1",
        phenotype="pheno1",
        gene_hits=hits,
        locus_pieces=pieces,
        locus_seqs=seqs,
        gene_states=np.array([GeneState.NORMAL.value, GeneState.NOVEL.value])
    )
    
    plotter = SerotypingResultPlotter()
    fig = plotter(res)
    print(fig.to_json()[:100])
    print("Success")

if __name__ == "__main__":
    main()
