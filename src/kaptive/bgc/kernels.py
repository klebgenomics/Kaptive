r"""Numba-accelerated numerical kernels for Biosynthetic Gene Cluster (BGC) feature extraction and sequence modeling.

High-performance parallelized Numba JIT functions for vectorizing ORF protein strobemer hashes,
applying sliding window contextual smoothing, calculating centroid emission likelihoods, and
computing optimal Viterbi state trajectories.

Functions:
    vectorize_orfs_kernel: Parallel kernel converting ORF strobemer hash records into L2-normalized feature vectors
        ([`vectorize_orfs_kernel`][kaptive.bgc.kernels.vectorize_orfs_kernel]).
    contextualize_features_kernel: Parallel kernel performing sliding-window spatial averaging of feature vectors
        ([`contextualize_features_kernel`][kaptive.bgc.kernels.contextualize_features_kernel]).
    archetype_emissions_kernel: Parallel kernel computing emission log-probabilities against background and archetype
        centroids ([`archetype_emissions_kernel`][kaptive.bgc.kernels.archetype_emissions_kernel]).
    architectural_viterbi_kernel: Dynamic programming Viterbi kernel finding maximum likelihood BGC architecture
        state paths ([`architectural_viterbi_kernel`][kaptive.bgc.kernels.architectural_viterbi_kernel]).
"""

import math

import numpy as np
from numba import njit, prange


@njit(parallel=True, fastmath=True)
def vectorize_orfs_kernel(records: np.ndarray, offsets: np.ndarray, num_seqs: int, feature_dim: int) -> np.ndarray:
    r"""Vectorize ORF strobemer hash records into L2-normalized feature vectors.

    Computes bag-of-strobemer count histograms for each sequence using modulo hashing,
    followed by L2 norm scaling across feature dimensions.

    Args:
        records (np.ndarray): 1D structured NumPy array containing strobemer hash entries with field `'hash'`.
        offsets (np.ndarray): 1D uint array of record offset boundaries per sequence.
        num_seqs (int): Total number of target sequences (ORFs) to process.
        feature_dim (int): Dimension of feature representation vector space.

    Returns:
        np.ndarray: 2D float32 array of shape `(num_seqs, feature_dim)` containing L2-normalized feature vectors.
    """
    features = np.zeros((num_seqs, feature_dim), dtype=np.float32)

    for i in prange(num_seqs):  # type: ignore
        start = offsets[i]
        end = offsets[i + 1] if i + 1 < num_seqs else len(records)

        for j in range(start, end):
            hash_val = records[j]["hash"]
            idx = hash_val % feature_dim
            features[i, idx] += 1.0

    # L2 normalize
    for i in prange(num_seqs):  # type: ignore
        norm = 0.0
        for j in range(feature_dim):
            norm += features[i, j] * features[i, j]
        if norm > 0:
            norm = math.sqrt(norm)
            for j in range(feature_dim):
                features[i, j] /= norm

    return features


@njit(parallel=True, fastmath=True)
def contextualize_features_kernel(features: np.ndarray, window_size: int) -> np.ndarray:
    r"""Apply sliding-window mean smoothing to ORF feature vectors.

    Aggregates neighboring feature representations within a genomic context window to incorporate
    spatial locus environment information.

    Args:
        features (np.ndarray): 2D float32 array of shape `(num_seqs, feature_dim)` containing input feature vectors.
        window_size (int): Half-width of sliding context window (number of adjacent ORFs to include on each side).

    Returns:
        np.ndarray: 2D float32 array of shape `(num_seqs, feature_dim)` containing context-smoothed feature vectors.
    """
    num_seqs, dim = features.shape
    out = np.zeros_like(features)

    for i in prange(num_seqs):  # type: ignore
        start = max(0, i - window_size)
        end = min(num_seqs, i + window_size + 1)
        count = end - start

        for j in range(start, end):
            for d in range(dim):
                out[i, d] += features[j, d]

        for d in range(dim):
            out[i, d] /= count

    return out


@njit(parallel=True, fastmath=True)
def archetype_emissions_kernel(
    features: np.ndarray, centroids: np.ndarray, bg_centroid: np.ndarray, novelty_score: float
) -> np.ndarray:
    r"""Calculate log emission probabilities for background, archetype, and novelty states.

    Evaluates negative squared Euclidean distance between contextualized ORF features and cluster
    archetype centroids, background profile, and a fixed novelty threshold.

    Args:
        features (np.ndarray): 2D float32 array of shape `(num_seqs, feature_dim)` containing context features.
        centroids (np.ndarray): 2D float32 array of shape `(num_archetypes, feature_dim)` storing archetype
            cluster centroids.
        bg_centroid (np.ndarray): 1D float32 array of length `feature_dim` representing non-BGC background
            feature centroid.
        novelty_score (float): Constant log-probability score assigned to unassigned/novel architectural state.

    Returns:
        np.ndarray: 2D float32 array of shape `(num_seqs, num_states)` storing log emission values,
            where `num_states = num_archetypes + 3`.
    """
    num_seqs, _ = features.shape
    num_archetypes = centroids.shape[0]
    num_states = num_archetypes + 3

    emissions = np.zeros((num_seqs, num_states), dtype=np.float32)

    for i in prange(num_seqs):  # type: ignore
        # 1. Background (PRE and POST)
        bg_dist = 0.0
        for d in range(features.shape[1]):
            diff = features[i, d] - bg_centroid[d]
            bg_dist += diff * diff

        bg_prob = -bg_dist
        emissions[i, 0] = bg_prob
        emissions[i, num_states - 1] = bg_prob

        # 2. Archetypes
        for k in range(num_archetypes):
            dist = 0.0
            for d in range(features.shape[1]):
                diff = features[i, d] - centroids[k, d]
                dist += diff * diff
            emissions[i, k + 1] = -dist

        # 3. NOVEL bridging state
        emissions[i, num_archetypes + 1] = novelty_score

    return emissions


@njit(fastmath=True)
def architectural_viterbi_kernel(emissions: np.ndarray, transitions: np.ndarray) -> tuple[np.ndarray, float]:
    r"""Compute optimal architectural state path using Viterbi dynamic programming.

    Solves maximum log-likelihood state trajectory over ORF emission probabilities and state
    transition matrix, subject to starting in PRE-BGC background state and terminating in
    POST-BGC background state.

    Args:
        emissions (np.ndarray): 2D float32 array of shape `(n_seqs, n_states)` containing log emission
            log-probabilities.
        transitions (np.ndarray): 2D float32 array of shape `(n_states, n_states)` containing log state
            transition probabilities.

    Returns:
        tuple[np.ndarray, float]: Tuple containing:
            - `path` (`np.ndarray`): 1D int32 array of length `n_seqs` containing optimal state index sequence.
            - `max_prob` (`float`): Total log path probability of optimal state trajectory, or `-inf` if `n_seqs == 0`.
    """
    n_seqs, n_states = emissions.shape

    if n_seqs == 0:
        return np.zeros(0, dtype=np.int32), float("-inf")

    viterbi = np.full((n_seqs, n_states), -np.inf, dtype=np.float32)
    backpointers = np.zeros((n_seqs, n_states), dtype=np.int32)

    # Init first step from PRE (State 0)
    for s in range(n_states):
        viterbi[0, s] = transitions[0, s] + emissions[0, s]

    for t in range(1, n_seqs):
        for s in range(n_states):
            max_tr_prob = -np.inf
            best_prev = -1

            for prev_s in range(n_states):
                prob = viterbi[t - 1, prev_s] + transitions[prev_s, s]
                if prob > max_tr_prob:
                    max_tr_prob = prob
                    best_prev = prev_s

            viterbi[t, s] = max_tr_prob + emissions[t, s]
            backpointers[t, s] = best_prev

    # Force termination in POST (State n_states - 1)
    max_prob = viterbi[-1, n_states - 1]
    best_last = n_states - 1

    path = np.zeros(n_seqs, dtype=np.int32)
    path[-1] = best_last

    for t in range(n_seqs - 2, -1, -1):
        path[t] = backpointers[t + 1, path[t + 1]]

    return path, max_prob
