r"""Biosynthetic Gene Cluster (BGC) locus predictor and state path container.

This module provides the Viterbi decoding engine
[`ArchitecturalPredictor`][kaptive.bgc.predictor.ArchitecturalPredictor]
and prediction output batch container [`BGCPredictions`][kaptive.bgc.predictor.BGCPredictions].

Classes:
    BGCPredictions: Container for predicted BGC locus boundaries, scores, and ORF state paths.
    ArchitecturalPredictor: Inference engine predicting BGC locus regions using Randstrobe features.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Self

import numpy as np
import numpy.typing as npt

from kaptive.bgc.annotate import AnnotationResult
from kaptive.bgc.kernels import (
    archetype_emissions_kernel,
    architectural_viterbi_kernel,
    contextualize_features_kernel,
    vectorize_orfs_kernel,
)
from kaptive.bgc.models import ArchitecturalModel
from kaptive.core.collections import BatchedContainer
from kaptive.core.kmers import RandstrobeIndex, _compute_query_offsets
from kaptive.core.seq import Sequences


@dataclass(slots=True, frozen=True)
class BGCPredictions(BatchedContainer[Any, "BGCPredictions"]):
    r"""Container for predicted BGC locus boundaries, scores, and state paths.

    Implements [`BatchedContainer`][kaptive.core.collections.BatchedContainer] interface
    for storing predictions across multiple contigs.

    Attributes:
        contig_names (npt.NDArray[np.object_]): Array of contig identifiers for each predicted locus.
        scores (npt.NDArray[np.float32]): Array of Viterbi path log-likelihood scores for each locus.
        orf_indices (list[npt.NDArray[np.uint32]]): List of arrays containing 0-based ORF indices
            belonging to each predicted BGC locus.
        paths (list[npt.NDArray[np.int32]]): List of arrays containing decoded HMM state sequences
            for predicted BGC ORFs.
        proteins_list (list[Sequences]): List of [`Sequences`][kaptive.core.seq.Sequences] containers
            holding translated amino acid sequences of predicted BGC ORFs.
    """

    contig_names: npt.NDArray[np.object_]
    scores: npt.NDArray[np.float32]
    orf_indices: list[npt.NDArray[np.uint32]]
    paths: list[npt.NDArray[np.int32]]
    proteins_list: list[Sequences]

    def __len__(self) -> int:
        r"""Return the number of predicted BGC loci in this batch.

        Returns:
            int: Number of predicted BGC locus entries.
        """
        return len(self.scores)

    def __getitem__(self, item: int | slice | npt.NDArray[Any] | list[int]) -> Any:  # noqa: ANN401
        r"""Index or slice predicted BGC locus entries.

        Args:
            item (int | slice | npt.NDArray | list): Integer index, slice object, boolean mask,
                or index array.

        Returns:
            tuple[str, float, npt.NDArray[np.uint32], npt.NDArray[np.int32], Sequences] | BGCPredictions:
                If `item` is an integer, returns a 5-element tuple `(contig_name, score, orf_indices, path, proteins)`.
                Otherwise, returns a sliced [`BGCPredictions`][kaptive.bgc.predictor.BGCPredictions] container.

        Raises:
            IndexError: If specified index is out of range.
        """
        if isinstance(item, int):
            return (
                self.contig_names[item],
                self.scores[item],
                self.orf_indices[item],
                self.paths[item],
                self.proteins_list[item],
            )
        return BGCPredictions(
            contig_names=self.contig_names[item],
            scores=self.scores[item],
            orf_indices=[self.orf_indices[i] for i in np.arange(len(self))[item]],
            paths=[self.paths[i] for i in np.arange(len(self))[item]],
            proteins_list=[self.proteins_list[i] for i in np.arange(len(self))[item]],
        )

    @classmethod
    def empty(cls) -> "BGCPredictions":
        r"""Construct an empty [`BGCPredictions`][kaptive.bgc.predictor.BGCPredictions] container.

        Returns:
            BGCPredictions: Empty container with zero-length arrays and empty lists.
        """
        return cls(
            contig_names=np.empty(0, dtype=object),
            scores=np.empty(0, dtype=np.float32),
            orf_indices=[],
            paths=[],
            proteins_list=[],
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenate multiple [`BGCPredictions`][kaptive.bgc.predictor.BGCPredictions] batches.

        Args:
            batches (list[BGCPredictions]): List of prediction containers to merge.

        Returns:
            BGCPredictions: Merged predictions container.
        """
        if not batches:
            return cls.empty()  # type: ignore
        return cls(
            contig_names=np.concatenate([b.contig_names for b in batches]),
            scores=np.concatenate([b.scores for b in batches]),
            orf_indices=[idx for b in batches for idx in b.orf_indices],
            paths=[p for b in batches for p in b.paths],
            proteins_list=[p for b in batches for p in b.proteins_list],
        )


class ArchitecturalPredictor:
    r"""Inference engine for predicting BGC locus regions from sequence annotations.

    Evaluates contig ORF translation features using Randstrobe indexing, feature smoothing,
    emission scoring against model centroids, and Viterbi state decoding using an
    [`ArchitecturalModel`][kaptive.bgc.models.ArchitecturalModel].

    Attributes:
        model (ArchitecturalModel): Model containing centroids and transition probabilities.
        window_size (int): Context window radius for feature vector smoothing.
        novelty_score (float): Penalty log-score assigned for unknown/unmatched feature profiles.
        dim (int): Feature vector dimensionality derived from model centroids shape.
        num_states (int): Number of HMM architectural states derived from model transition matrix shape.
    """

    def __init__(
        self,
        model: ArchitecturalModel,
        window_size: int = 1,
        novelty_score: float = -1.0,
    ) -> None:
        r"""Initialize the architectural predictor engine.

        Args:
            model (ArchitecturalModel): Trained architectural model instance.
            window_size (int): Context window radius for smoothing ORF feature profiles.
                Defaults to `1`.
            novelty_score (float): Penalty score for unmatched feature profiles.
                Defaults to `-1.0`.
        """
        self.model = model
        self.window_size = window_size
        self.novelty_score = novelty_score
        self.dim = self.model.centroids.shape[1]
        self.num_states = self.model.transitions.shape[0]

    def predict_contig(self, proteins: Sequences) -> tuple[npt.NDArray[np.int32], float]:
        r"""Predict BGC architectural state path and score for a single contig's protein sequences.

        Args:
            proteins (Sequences): Translated protein sequences for all predicted ORFs in a contig.

        Returns:
            tuple[npt.NDArray[np.int32], float]: A 2-element tuple `(path, score)`:
                - `path`: Array of 0-based architectural state assignments for each ORF.
                - `score`: Total Viterbi log-likelihood score for the predicted optimal path.
        """
        num_orfs = len(proteins)
        if num_orfs == 0:
            return np.zeros(0, dtype=np.int32), float("-inf")

        idx = RandstrobeIndex.build(proteins, sort_by_hash=False)
        if len(idx) == 0:
            return np.zeros(num_orfs, dtype=np.int32), float("-inf")

        offsets = _compute_query_offsets(idx.records, idx.n_seqs)

        X_raw = vectorize_orfs_kernel(idx.records, offsets, num_orfs, self.dim)
        X_smooth = contextualize_features_kernel(X_raw, self.window_size)
        emissions = archetype_emissions_kernel(
            X_smooth, self.model.centroids, self.model.bg_centroid, self.novelty_score
        )

        # Mask emissions for genes at the edge of the contig (they may be truncated)
        if num_orfs > 0:
            emissions[0, :] = 0.0
            emissions[-1, :] = 0.0

        path, score = architectural_viterbi_kernel(emissions, self.model.transitions)
        return path, score

    def predict(self, annotation: AnnotationResult) -> BGCPredictions:
        r"""Predict BGC locus regions across all contigs in an annotation result.

        Args:
            annotation (AnnotationResult): Annotation output containing contig names and predicted genes
                ([`AnnotationResult`][kaptive.bgc.annotate.AnnotationResult]).

        Returns:
            BGCPredictions: Container holding predicted BGC locus boundaries, scores, ORF indices,
                and state paths ([`BGCPredictions`][kaptive.bgc.predictor.BGCPredictions]).
        """
        contig_names, scores, orf_indices, paths, proteins_list = [], [], [], [], []

        for c_idx, c_name in enumerate(annotation.contig_names):
            mask = annotation.genes.contig_indices == c_idx
            if not np.any(mask):
                continue

            contig_genes = annotation.genes[mask]
            path, score = self.predict_contig(contig_genes.translations)

            if score > -np.inf:
                # Identify BGC ORFs (states 1 to num_states - 2)
                bgc_mask = (path > 0) & (path < self.num_states - 1)
                bgc_inds = np.where(bgc_mask)[0]

                if len(bgc_inds) > 0:
                    contig_names.append(c_name)
                    scores.append(score)
                    orf_indices.append(bgc_inds.astype(np.uint32))
                    paths.append(path[bgc_inds])
                    proteins_list.append(contig_genes.translations[bgc_inds])

        return BGCPredictions(
            contig_names=np.array(contig_names, dtype=object),
            scores=np.array(scores, dtype=np.float32),
            orf_indices=orf_indices,
            paths=paths,
            proteins_list=proteins_list,
        )
