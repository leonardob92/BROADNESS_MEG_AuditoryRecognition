"""Omnibus repeated-measures temporal cluster permutation test."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import f as f_distribution

from ._temporal_clusters import (
    TemporalClusterSettings,
    TemporalClusterStats,
    _cluster_results,
    _clusters,
    _options,
    _repeated_f,
    _series,
)


def BROADNESS_TemporalClusterCorrection_ANOVA(
    network_time_series: ArrayLike,
    alpha: float = 0.05,
    *,
    permutations: int = 5000,
    clusteralpha: float = 0.05,
    randomseed: int = 1,
    time: ArrayLike | None = None,
) -> tuple[NDArray[np.bool_], NDArray[np.float64], float, TemporalClusterStats]:
    """Test the condition effect with whole-series label permutations.

    A different condition-label permutation is drawn for each participant
    and applied unchanged over time. Each network has a separate maximum
    positive-F cluster null distribution.
    """
    values = _series(network_time_series)
    n_time, n_networks, n_conditions, n_participants = values.shape
    if n_time < 2 or n_networks < 1 or n_conditions < 2 or n_participants < 2:
        raise ValueError("need at least 2 time points, 1 network, 2 conditions and 2 participants")
    level, forming, n_permutations, seed, times = _options(
        alpha, clusteralpha, permutations, randomseed, time, n_time
    )
    df_condition = n_conditions - 1
    df_error = (n_participants - 1) * df_condition
    threshold = float(f_distribution.ppf(1 - forming, df_condition, df_error))
    shape = (n_time, n_networks)
    significant = np.zeros(shape, dtype=bool)
    corrected = np.full(shape, np.nan)
    f_statistics = np.zeros(shape)
    cluster_labels = np.zeros(shape, dtype=np.int64)
    cluster_tables = np.empty((n_networks,), dtype=object)
    null_maxima = np.zeros((n_permutations, n_networks))
    generator = np.random.default_rng(seed)

    for network in range(n_networks):
        observed = values[:, network, :, :]
        observed_f = _repeated_f(observed)
        labels, masses = _clusters(observed_f, threshold, two_sided=False)
        f_statistics[:, network] = observed_f
        cluster_labels[:, network] = labels

        for permutation in range(n_permutations):
            permuted = np.empty_like(observed)
            for participant in range(n_participants):
                order = generator.permutation(n_conditions)
                permuted[:, :, participant] = observed[:, order, participant]
            permuted_f = _repeated_f(permuted)
            _, permuted_masses = _clusters(permuted_f, threshold, two_sided=False)
            if permuted_masses.size:
                null_maxima[permutation, network] = np.max(permuted_masses)

        mask, p_map, table = _cluster_results(
            observed_f,
            labels,
            masses,
            null_maxima[:, network],
            level,
            times,
            f_statistic=True,
        )
        significant[:, network] = mask
        corrected[:, network] = p_map
        cluster_tables[network] = table

    settings = TemporalClusterSettings(
        Alpha=level,
        ClusterFormingAlpha=forming,
        Permutations=n_permutations,
        RandomSeed=seed,
        Time=times,
        CorrectionFamily="All time points, separately per network (omnibus condition effect)",
        Resampling="Complete condition-label permutations within participants",
        InferenceCaveat="Condition exchangeability required; cluster-level inference only; no across-network correction",
    )
    details = TemporalClusterStats(
        Settings=settings,
        FStatistics=f_statistics,
        ClusterFormingThreshold=threshold,
        DegreesFreedom=np.array([df_condition, df_error], dtype=float),
        ClusterLabels=cluster_labels,
        Clusters=cluster_tables,
        NullMaximumClusterMass=null_maxima,
        ClusterCorrectedP=corrected,
        Significant=significant,
    )
    return significant, corrected, level, details


__all__ = ["BROADNESS_TemporalClusterCorrection_ANOVA"]
