"""Independent-group Welch ANOVA temporal cluster wild bootstrap."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import f as f_distribution

from ._temporal_clusters import (
    TemporalClusterSettings,
    TemporalClusterStats,
    _cluster_results,
    _clusters,
    _integer_vector,
    _options,
    _series,
    _welch_f,
)


def BROADNESS_TemporalClusterCorrection_GroupANOVA(
    network_time_series: ArrayLike,
    group_labels: ArrayLike,
    alpha: float = 0.05,
    *,
    permutations: int = 5000,
    clusteralpha: float = 0.05,
    randomseed: int = 1,
    time: ArrayLike | None = None,
) -> tuple[NDArray[np.bool_], NDArray[np.float64], float, TemporalClusterStats]:
    """Test three or more independent groups per network and condition.

    The null uses group-centered participant residuals multiplied by one
    Rademacher sign per participant over the entire time series. This is an
    approximate heteroscedasticity-aware wild bootstrap.
    """
    values = _series(network_time_series)
    n_time, n_networks, n_conditions, n_participants = values.shape
    if n_time < 2 or n_networks < 1 or n_conditions < 1:
        raise ValueError("need at least 2 time points, 1 network and 1 condition")
    labels = _integer_vector(group_labels, "group_labels", n_participants)
    group_ids = np.unique(labels)
    if group_ids.size < 3:
        raise ValueError("overall Welch ANOVA requires at least three groups")
    group_indices = [np.flatnonzero(labels == group) for group in group_ids]
    if any(indices.size < 2 for indices in group_indices):
        raise ValueError("every group needs at least two participants")
    level, forming, n_permutations, seed, times = _options(
        alpha, clusteralpha, permutations, randomseed, time, n_time
    )
    shape = (n_time, n_networks, n_conditions)
    significant = np.zeros(shape, dtype=bool)
    corrected = np.full(shape, np.nan)
    f_statistics = np.zeros(shape)
    cluster_labels = np.zeros(shape, dtype=np.int64)
    cluster_threshold = np.full(shape, np.nan)
    cluster_tables = np.empty((n_networks, n_conditions), dtype=object)
    null_maxima = np.zeros((n_permutations, n_networks, n_conditions))
    generator = np.random.default_rng(seed)
    df1 = group_ids.size - 1

    for network in range(n_networks):
        for condition in range(n_conditions):
            groups = [values[:, network, condition, indices] for indices in group_indices]
            residuals = [group - np.mean(group, axis=1, keepdims=True) for group in groups]
            observed_f, observed_df2 = _welch_f(groups)
            threshold = f_distribution.ppf(1 - forming, df1, observed_df2)
            f_statistics[:, network, condition] = observed_f
            cluster_threshold[:, network, condition] = threshold
            observed_labels, masses = _clusters(observed_f, threshold, two_sided=False)
            cluster_labels[:, network, condition] = observed_labels

            for permutation in range(n_permutations):
                boot_groups = [
                    residual * generator.choice((-1.0, 1.0), size=residual.shape[1])
                    for residual in residuals
                ]
                boot_f, boot_df2 = _welch_f(boot_groups)
                boot_threshold = f_distribution.ppf(1 - forming, df1, boot_df2)
                _, boot_masses = _clusters(boot_f, boot_threshold, two_sided=False)
                if boot_masses.size:
                    null_maxima[permutation, network, condition] = np.max(boot_masses)

            mask, p_map, table = _cluster_results(
                observed_f,
                observed_labels,
                masses,
                null_maxima[:, network, condition],
                level,
                times,
                f_statistic=True,
            )
            significant[:, network, condition] = mask
            corrected[:, network, condition] = p_map
            cluster_tables[network, condition] = table

    settings = TemporalClusterSettings(
        Alpha=level,
        ClusterFormingAlpha=forming,
        Permutations=n_permutations,
        RandomSeed=seed,
        Time=times,
        GroupLabels=labels,
        GroupIDs=group_ids,
        CorrectionFamily="All time points, separately per network and condition",
        Resampling="Within-group centered-residual Rademacher wild bootstrap",
        InferenceCaveat="Approximate cluster-level inference; validate with null simulations, especially for small or unequal groups",
    )
    details = TemporalClusterStats(
        Settings=settings,
        FStatistics=f_statistics,
        ClusterFormingThreshold=cluster_threshold,
        ClusterLabels=cluster_labels,
        Clusters=cluster_tables,
        NullMaximumClusterMass=null_maxima,
        ClusterCorrectedP=corrected,
        Significant=significant,
    )
    return significant, corrected, level, details


__all__ = ["BROADNESS_TemporalClusterCorrection_GroupANOVA"]
