r"""Training engine for Biosynthetic Gene Cluster (BGC) architectural models.

This module provides the [`ArchitecturalTrainer`][kaptive.bgc.trainer.ArchitecturalTrainer]
class, which trains hidden Markov model (HMM) style architectural models
([`ArchitecturalModel`][kaptive.bgc.models.ArchitecturalModel]) from training
genomic loci ([`LocusData`][kaptive.compare.LocusData]) and background protein sequences.

The training pipeline:
1. Concatenates all core locus proteins and extracts randstrobe feature vectors using
   [`vectorize_orfs_kernel`][kaptive.bgc.kernels.vectorize_orfs_kernel].
2. Clusters core protein features into a discrete set of archetype centroids using
   mini-batch $k$-means (`sklearn.cluster.MiniBatchKMeans`).
3. Computes background protein feature centroid for reference non-BGC background emission scoring.
4. Learns state transition frequencies between background (`PRE`), archetypes ($1 \dots K$),
   `NOVEL` gene states, and termination (`POST`).
5. Applies Laplacian smoothing and calculates transition log-probabilities.

Classes:
    [`ArchitecturalTrainer`][kaptive.bgc.trainer.ArchitecturalTrainer]: Trainer engine for
        building architectural BGC prediction models.
"""

from collections.abc import Sequence

import numpy as np
from sklearn.cluster import MiniBatchKMeans

from kaptive.bgc.kernels import vectorize_orfs_kernel
from kaptive.bgc.models import ArchitecturalModel
from kaptive.compare import LocusData
from kaptive.core.kmers import RandstrobeIndex, _compute_query_offsets
from kaptive.core.seq import Sequences


class ArchitecturalTrainer:
    r"""Trainer engine for building BGC hidden Markov model architectures.

    Learns protein archetype centroids and state transition log-probabilities
    from known BGC loci ([`LocusData`][kaptive.compare.LocusData]) and background
    proteins ([`Sequences`][kaptive.core.seq.Sequences]). The resulting model can be
    saved and used by [`ArchitecturalPredictor`][kaptive.bgc.predictor.ArchitecturalPredictor]
    for Viterbi decoding of novel gene clusters.

    Attributes:
        num_archetypes (int): Number of protein archetype clusters to fit via $k$-means.
        feature_dim (int): Dimensionality of the randstrobe hash feature space.
        partial_edge_tolerance (float): Laplacian smoothing weight added to entry/exit
            transitions for handling partial or fragmented contig edges.
    """

    def __init__(
        self,
        num_archetypes: int = 25,
        feature_dim: int = 1024,
        partial_edge_tolerance: float = 5.0,
    ) -> None:
        r"""Initialize the architectural trainer with clustering and HMM parameters.

        Args:
            num_archetypes (int): Number of protein archetype clusters ($K$) to fit. Defaults to `25`.
            feature_dim (int): Dimension of the randstrobe feature vector space ($D$). Defaults to `1024`.
            partial_edge_tolerance (float): Smoothing weight applied to transitions from `PRE` to archetypes
                and archetypes to `POST` to accommodate incomplete contig boundaries. Defaults to `5.0`.
        """
        self.num_archetypes = num_archetypes
        self.feature_dim = feature_dim
        self.partial_edge_tolerance = partial_edge_tolerance

    def train(
        self,
        loci: Sequence[LocusData],
        background_proteins: Sequences | None = None,
    ) -> ArchitecturalModel:
        r"""Train an architectural model from annotated loci and optional background proteins.

        Extracts features for all protein sequences in `loci`, clusters them into `num_archetypes`
        centroids using `MiniBatchKMeans`, builds a background centroid from `background_proteins`,
        and constructs a transition log-probability matrix with structural constraints and Laplacian smoothing.

        Args:
            loci (Sequence[LocusData]): Sequence of annotated BGC loci containing protein sequences
                in [`LocusData`][kaptive.compare.LocusData].
            background_proteins (Sequences | None): Optional background protein sequences used to derive
                non-BGC background emission profile. If `None` or empty, defaults to a zero vector.
                Defaults to `None`.

        Returns:
            ArchitecturalModel: The trained model instance containing centroids, background centroid,
                and log-transition matrix.

        Raises:
            ValueError: If `loci` is empty or contains fewer protein sequences than `num_archetypes`.
        """
        # 1. Gather all core locus sequences
        all_seqs = Sequences.concat([locus.proteins for locus in loci])

        # 2. Vectorize core genes
        idx = RandstrobeIndex.build(all_seqs, sort_by_hash=False)
        offsets = _compute_query_offsets(idx.records, idx.n_seqs)
        X_core = vectorize_orfs_kernel(idx.records, offsets, len(all_seqs), self.feature_dim)

        # 3. Vectorize background genes
        if background_proteins is not None and len(background_proteins) > 0:
            bg_idx = RandstrobeIndex.build(background_proteins, sort_by_hash=False)
            bg_offsets = _compute_query_offsets(bg_idx.records, bg_idx.n_seqs)
            X_bg = vectorize_orfs_kernel(bg_idx.records, bg_offsets, len(background_proteins), self.feature_dim)
            bg_centroid = X_bg.mean(axis=0)
        else:
            bg_centroid = np.zeros(self.feature_dim, dtype=np.float32)

        # 4. Cluster core genes
        kmeans = MiniBatchKMeans(n_clusters=self.num_archetypes, random_state=42, n_init="auto")
        archetypes = kmeans.fit_predict(X_core)
        centroids = kmeans.cluster_centers_

        # Map global sequences back to archetypes (shifting by 1 because 0 is PRE)
        gene_state_map = archetypes + 1

        # 5. Learn Transition Matrix
        num_states = self.num_archetypes + 3
        transitions = np.zeros((num_states, num_states), dtype=np.float32)

        PRE, NOVEL, POST = 0, self.num_archetypes + 1, self.num_archetypes + 2

        offset = 0
        for locus in loci:
            length = len(locus.proteins)
            if length == 0:
                continue

            prev_state = PRE
            for i in range(length):
                curr_state = gene_state_map[offset + i]
                transitions[prev_state, curr_state] += 1.0
                prev_state = curr_state

            transitions[prev_state, POST] += 1.0
            offset += length

        # Apply Laplacian smoothing to avoid -inf transitions
        transitions += 0.1

        # Explicitly weight background self-loops and NOVEL bridges
        transitions[PRE, PRE] += 100.0
        transitions[POST, POST] += 100.0
        transitions[1:NOVEL, NOVEL] += 5.0
        transitions[NOVEL, 1:NOVEL] += 5.0
        transitions[NOVEL, NOVEL] += 10.0

        # Allow entry/exit anywhere if contig is fragmented
        transitions[PRE, 1:NOVEL] += self.partial_edge_tolerance
        transitions[1:NOVEL, POST] += self.partial_edge_tolerance

        # Convert to log probabilities
        row_sums = transitions.sum(axis=1, keepdims=True)
        trans_log_probs = np.log(transitions / row_sums)

        # Hard structural constraints (Cannot escape POST)
        trans_log_probs[POST, :] = -np.inf
        trans_log_probs[POST, POST] = 0.0

        return ArchitecturalModel(
            centroids=np.ascontiguousarray(centroids.astype(np.float32)),
            bg_centroid=np.ascontiguousarray(bg_centroid.astype(np.float32)),
            transitions=np.ascontiguousarray(trans_log_probs.astype(np.float32)),
        )
