"""Pairwise independent-group Welch temporal cluster wild bootstrap."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import t as t_distribution

from ._temporal_clusters import (
    TemporalClusterSettings,
    TemporalClusterStats,
    _cluster_results,
    _clusters,
    _integer_vector,
    _options,
    _pairs,
    _series,
    _welch_t,
)


def BROADNESS_TemporalClusterCorrection_Groups(
    network_time_series: ArrayLike,
    group_labels: ArrayLike,
    group_pairs: ArrayLike,
    alpha: float = 0.05,
    *,
    permutations: int = 5000,
    clusteralpha: float = 0.05,
    randomseed: int = 1,
    time: ArrayLike | None = None,
) -> tuple[NDArray[np.bool_], NDArray[np.float64], float, TemporalClusterStats]:
    """Compare independent groups separately per network, condition, and pair.

    The null uses group-centered participant residuals multiplied by one
    Rademacher sign per participant over the entire time series. This is an
    approximate wild bootstrap, not an exact label permutation.
    """
    values = _series(network_time_series)
    n_time, n_networks, n_conditions, n_participants = values.shape
    if n_time < 2 or n_networks < 1 or n_conditions < 1 or n_participants < 4:
        raise ValueError("need at least 2 time points, 1 network, 1 condition and 4 participants")
    labels = _integer_vector(group_labels, "group_labels", n_participants)
    pairs = _pairs(group_pairs, "group_pairs", np.unique(labels))
    group_indices = []
    for first, second in pairs:
        first_indices = np.flatnonzero(labels == first)
        second_indices = np.flatnonzero(labels == second)
        if first_indices.size < 2 or second_indices.size < 2:
            raise ValueError("each requested group needs at least two participants")
        group_indices.append((first_indices, second_indices))
    level, forming, n_permutations, seed, times = _options(
        alpha, clusteralpha, permutations, randomseed, time, n_time
    )
    n_pairs = pairs.shape[0]
    shape = (n_time, n_networks, n_conditions, n_pairs)
    significant = np.zeros(shape, dtype=bool)
    corrected = np.full(shape, np.nan)
    t_statistics = np.zeros(shape)
    cluster_labels = np.zeros(shape, dtype=np.int64)
    cluster_threshold = np.full(shape, np.nan)
    cluster_tables = np.empty((n_networks, n_conditions, n_pairs), dtype=object)
    null_maxima = np.zeros((n_permutations, n_networks, n_conditions, n_pairs))
    generator = np.random.default_rng(seed)

    for network in range(n_networks):
        for condition in range(n_conditions):
            for pair_index, (first_indices, second_indices) in enumerate(group_indices):
                x = values[:, network, condition, first_indices]
                y = values[:, network, condition, second_indices]
                observed_t, observed_df = _welch_t(x, y)
                threshold = t_distribution.ppf(1 - forming / 2, observed_df)
                t_statistics[:, network, condition, pair_index] = observed_t
                cluster_threshold[:, network, condition, pair_index] = threshold
                observed_labels, masses = _clusters(observed_t, threshold, two_sided=True)
                cluster_labels[:, network, condition, pair_index] = observed_labels

                residual_x = x - np.mean(x, axis=1, keepdims=True)
                residual_y = y - np.mean(y, axis=1, keepdims=True)
                for permutation in range(n_permutations):
                    signs_x = generator.choice((-1.0, 1.0), size=x.shape[1])
                    signs_y = generator.choice((-1.0, 1.0), size=y.shape[1])
                    boot_t, boot_df = _welch_t(
                        residual_x * signs_x, residual_y * signs_y
                    )
                    boot_threshold = t_distribution.ppf(1 - forming / 2, boot_df)
                    _, boot_masses = _clusters(boot_t, boot_threshold, two_sided=True)
                    if boot_masses.size:
                        null_maxima[permutation, network, condition, pair_index] = np.max(boot_masses)

                mask, p_map, table = _cluster_results(
                    observed_t,
                    observed_labels,
                    masses,
                    null_maxima[:, network, condition, pair_index],
                    level,
                    times,
                    f_statistic=False,
                )
                significant[:, network, condition, pair_index] = mask
                corrected[:, network, condition, pair_index] = p_map
                cluster_tables[network, condition, pair_index] = table

    settings = TemporalClusterSettings(
        Alpha=level,
        ClusterFormingAlpha=forming,
        Permutations=n_permutations,
        RandomSeed=seed,
        Time=times,
        GroupPairs=pairs,
        GroupLabels=labels,
        CorrectionFamily="All time points, separately per network, condition and group pair",
        Resampling="Within-group centered-residual Rademacher wild bootstrap",
        InferenceCaveat="Approximate cluster-level inference; validate with null simulations, especially for small or unequal groups",
    )
    details = TemporalClusterStats(
        Settings=settings,
        TStatistics=t_statistics,
        ClusterFormingThreshold=cluster_threshold,
        ClusterLabels=cluster_labels,
        Clusters=cluster_tables,
        NullMaximumClusterMass=null_maxima,
        ClusterCorrectedP=corrected,
        Significant=significant,
    )
    return significant, corrected, level, details


__all__ = ["BROADNESS_TemporalClusterCorrection_Groups"]
