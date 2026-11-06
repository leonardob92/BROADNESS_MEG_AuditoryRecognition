"""Paired temporal cluster-mass permutation test for condition contrasts."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import t as t_distribution

from ._temporal_clusters import (
    TemporalClusterSettings,
    TemporalClusterStats,
    _cluster_results,
    _clusters,
    _options,
    _paired_t,
    _pairs,
    _series,
)


def BROADNESS_TemporalClusterCorrection(
    network_time_series: ArrayLike,
    condition_pairs: ArrayLike,
    alpha: float = 0.05,
    *,
    permutations: int = 5000,
    clusteralpha: float = 0.05,
    randomseed: int = 1,
    time: ArrayLike | None = None,
) -> tuple[NDArray[np.bool_], NDArray[np.float64], float, TemporalClusterStats]:
    """Test paired conditions with complete-series participant sign flips.

    Input axes are time × selected networks × conditions × participants.
    Condition pairs use one-based numbers. Each network and pair has its own
    maximum-cluster null distribution across all time samples. Corrected p
    values belong to clusters and are repeated over their member samples.
    """
    values = _series(network_time_series)
    n_time, n_networks, n_conditions, n_participants = values.shape
    if n_time < 2 or n_networks < 1 or n_conditions < 2 or n_participants < 2:
        raise ValueError("need at least 2 time points, 1 network, 2 conditions and 2 participants")
    pairs = _pairs(condition_pairs, "condition_pairs", np.arange(1, n_conditions + 1))
    level, forming, n_permutations, seed, times = _options(
        alpha, clusteralpha, permutations, randomseed, time, n_time
    )
    n_pairs = pairs.shape[0]
    shape = (n_time, n_networks, n_pairs)
    significant = np.zeros(shape, dtype=bool)
    corrected = np.full(shape, np.nan)
    t_statistics = np.zeros(shape)
    cluster_labels = np.zeros(shape, dtype=np.int64)
    cluster_tables = np.empty((n_networks, n_pairs), dtype=object)
    null_maxima = np.zeros((n_permutations, n_networks, n_pairs))
    threshold = float(t_distribution.ppf(1 - forming / 2, n_participants - 1))
    generator = np.random.default_rng(seed)

    for network in range(n_networks):
        for pair_index, (first, second) in enumerate(pairs):
            difference = values[:, network, first - 1, :] - values[:, network, second - 1, :]
            observed_t = _paired_t(difference)
            labels, masses = _clusters(observed_t, threshold, two_sided=True)
            t_statistics[:, network, pair_index] = observed_t
            cluster_labels[:, network, pair_index] = labels

            for permutation in range(n_permutations):
                signs = generator.choice((-1.0, 1.0), size=n_participants)
                permuted_t = _paired_t(difference * signs)
                _, permuted_masses = _clusters(permuted_t, threshold, two_sided=True)
                if permuted_masses.size:
                    null_maxima[permutation, network, pair_index] = np.max(permuted_masses)

            mask, p_map, table = _cluster_results(
                observed_t,
                labels,
                masses,
                null_maxima[:, network, pair_index],
                level,
                times,
                f_statistic=False,
            )
            significant[:, network, pair_index] = mask
            corrected[:, network, pair_index] = p_map
            cluster_tables[network, pair_index] = table

    settings = TemporalClusterSettings(
        Alpha=level,
        ClusterFormingAlpha=forming,
        Permutations=n_permutations,
        RandomSeed=seed,
        Time=times,
        ConditionPairs=pairs,
        CorrectionFamily="All time points, separately per network and condition pair",
        Resampling="Paired-difference participant sign flips",
        InferenceCaveat="Sign-symmetric paired differences required; no correction across pairs or networks",
    )
    details = TemporalClusterStats(
        Settings=settings,
        TStatistics=t_statistics,
        ClusterFormingThreshold=threshold,
        ClusterLabels=cluster_labels,
        Clusters=cluster_tables,
        NullMaximumClusterMass=null_maxima,
        ClusterCorrectedP=corrected,
        Significant=significant,
    )
    return significant, corrected, level, details


__all__ = ["BROADNESS_TemporalClusterCorrection"]
