"""Shared validation, statistics, and cluster packaging for temporal tests."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
import pandas as pd

from .BROADNESS_FDRCorrection import _probability


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True, slots=True)
class TemporalClusterSettings:
    Alpha: float
    ClusterFormingAlpha: float
    Permutations: int
    RandomSeed: int
    Time: FloatArray | None
    CorrectionFamily: str
    Resampling: str
    InferenceCaveat: str
    ConditionPairs: IntArray | None = None
    GroupPairs: IntArray | None = None
    GroupLabels: IntArray | None = None
    GroupIDs: IntArray | None = None


@dataclass(slots=True)
class TemporalClusterStats:
    """MATLAB-style temporal cluster details; ``Clusters`` holds DataFrames."""

    Settings: TemporalClusterSettings
    ClusterFormingThreshold: float | FloatArray
    ClusterLabels: IntArray
    Clusters: NDArray[np.object_]
    NullMaximumClusterMass: FloatArray
    ClusterCorrectedP: FloatArray
    Significant: NDArray[np.bool_]
    TStatistics: FloatArray | None = None
    FStatistics: FloatArray | None = None
    DegreesFreedom: FloatArray | None = None


def _series(data: ArrayLike) -> FloatArray:
    values = np.asarray(data)
    if values.ndim != 4 or not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise ValueError("network_time_series must be a real time x networks x conditions x participants array")
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("network_time_series must be finite")
    return values


def _options(
    alpha: float,
    clusteralpha: float,
    permutations: int,
    randomseed: int,
    time: ArrayLike | None,
    n_time: int,
) -> tuple[float, float, int, int, FloatArray | None]:
    level = _probability(alpha, "alpha")
    cluster_level = _probability(clusteralpha, "clusteralpha")
    if isinstance(permutations, (bool, np.bool_)) or not isinstance(
        permutations, (int, np.integer)
    ) or permutations < 1:
        raise ValueError("permutations must be a positive integer")
    if isinstance(randomseed, (bool, np.bool_)) or not isinstance(
        randomseed, (int, np.integer)
    ) or randomseed < 0:
        raise ValueError("randomseed must be a nonnegative integer")
    if time is None:
        time_values = None
    else:
        raw = np.asarray(time)
        if (
            raw.ndim not in (1, 2)
            or (raw.ndim == 2 and 1 not in raw.shape)
            or raw.size != n_time
            or not np.issubdtype(raw.dtype, np.number)
            or np.iscomplexobj(raw)
        ):
            raise ValueError("time must be a real vector with one value per sample")
        time_values = np.asarray(raw, dtype=float).reshape(-1)
        if not np.all(np.isfinite(time_values)) or np.any(np.diff(time_values) <= 0):
            raise ValueError("time must be finite and strictly increasing")
    return level, cluster_level, int(permutations), int(randomseed), time_values


def _integer_vector(values: ArrayLike, name: str, expected: int) -> IntArray:
    raw = np.asarray(values)
    if (
        raw.ndim not in (1, 2)
        or (raw.ndim == 2 and 1 not in raw.shape)
        or raw.size != expected
        or not np.issubdtype(raw.dtype, np.number)
        or np.iscomplexobj(raw)
        or not np.all(np.isfinite(raw))
        or np.any(raw < 1)
        or np.any(raw != np.floor(raw))
    ):
        raise ValueError(f"{name} must contain one positive integer per participant")
    return np.asarray(raw, dtype=np.int64).reshape(-1)


def _pairs(values: ArrayLike, name: str, allowed: IntArray) -> IntArray:
    raw = np.asarray(values)
    if (
        raw.ndim != 2
        or raw.shape[1] != 2
        or raw.shape[0] < 1
        or not np.issubdtype(raw.dtype, np.number)
        or np.iscomplexobj(raw)
        or not np.all(np.isfinite(raw))
        or np.any(raw != np.floor(raw))
    ):
        raise ValueError(f"{name} must contain valid distinct nonduplicated pairs")
    pairs = np.asarray(raw, dtype=np.int64)
    canonical = np.sort(pairs, axis=1)
    if (
        np.any(pairs[:, 0] == pairs[:, 1])
        or not np.all(np.isin(pairs, allowed))
        or np.unique(canonical, axis=0).shape[0] != pairs.shape[0]
    ):
        raise ValueError(f"{name} must contain valid distinct nonduplicated pairs")
    return pairs


def _clusters(
    statistics: FloatArray, threshold: float | FloatArray, *, two_sided: bool
) -> tuple[IntArray, FloatArray]:
    """Label runs, positive first then negative as in the MATLAB functions."""
    statistic = np.asarray(statistics, dtype=float).reshape(-1)
    limit = np.broadcast_to(np.asarray(threshold, dtype=float), statistic.shape)
    labels = np.zeros(statistic.size, dtype=np.int64)
    masses: list[float] = []
    for direction in ((1, -1) if two_sided else (1,)):
        active = np.isfinite(limit) & (direction * statistic >= limit)
        changes = np.diff(np.concatenate(([False], active, [False])).astype(int))
        starts = np.flatnonzero(changes == 1)
        ends = np.flatnonzero(changes == -1)
        for start, stop in zip(starts, ends):
            labels[start:stop] = len(masses) + 1
            masses.append(float(np.sum(np.abs(statistic[start:stop])) if two_sided else np.sum(statistic[start:stop])))
    return labels, np.asarray(masses, dtype=float)


def _paired_t(differences: FloatArray) -> FloatArray:
    n_participants = differences.shape[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.mean(differences, axis=1) / (
            np.std(differences, axis=1, ddof=1) / np.sqrt(n_participants)
        )
    return np.where(np.isnan(result), 0.0, result)


def _repeated_f(data: FloatArray) -> FloatArray:
    """Classical one-factor repeated-measures F at every sample."""
    _, n_conditions, n_participants = data.shape
    grand = np.mean(data, axis=(1, 2), keepdims=True)
    condition_means = np.mean(data, axis=2, keepdims=True)
    participant_means = np.mean(data, axis=1, keepdims=True)
    ss_condition = n_participants * np.sum((condition_means - grand) ** 2, axis=(1, 2))
    ss_participant = n_conditions * np.sum((participant_means - grand) ** 2, axis=(1, 2))
    ss_total = np.sum((data - grand) ** 2, axis=(1, 2))
    ss_error = np.maximum(0.0, ss_total - ss_condition - ss_participant)
    df1 = n_conditions - 1
    df2 = (n_participants - 1) * df1
    with np.errstate(divide="ignore", invalid="ignore"):
        result = (ss_condition / df1) / (ss_error / df2)
    return np.where(np.isnan(result), 0.0, result)


def _welch_t(x: FloatArray, y: FloatArray) -> tuple[FloatArray, FloatArray]:
    n1, n2 = x.shape[1], y.shape[1]
    v1 = np.var(x, axis=1, ddof=1) / n1
    v2 = np.var(y, axis=1, ddof=1) / n2
    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = (np.mean(x, axis=1) - np.mean(y, axis=1)) / np.sqrt(v1 + v2)
        df = (v1 + v2) ** 2 / (v1**2 / (n1 - 1) + v2**2 / (n2 - 1))
    statistic = np.where(np.isfinite(statistic), statistic, 0.0)
    df = np.where(np.isfinite(df) & (df > 0), df, n1 + n2 - 2)
    return statistic, df


def _welch_f(groups: list[FloatArray]) -> tuple[FloatArray, FloatArray]:
    n_time = groups[0].shape[0]
    sizes = np.asarray([group.shape[1] for group in groups], dtype=float)
    means = np.column_stack([np.mean(group, axis=1) for group in groups])
    variances = np.column_stack([np.var(group, axis=1, ddof=1) for group in groups])
    valid = np.all(np.isfinite(variances) & (variances > 0), axis=1)
    statistic = np.zeros(n_time)
    df2 = np.full(n_time, np.nan)
    if not np.any(valid):
        return statistic, df2
    k = len(groups)
    weights = sizes / variances[valid]
    weight_total = np.sum(weights, axis=1)
    weighted_mean = np.sum(weights * means[valid], axis=1) / weight_total
    correction = np.sum(
        (1 - weights / weight_total[:, None]) ** 2 / (sizes - 1), axis=1
    )
    df2[valid] = (k**2 - 1) / (3 * correction)
    numerator = np.sum(weights * (means[valid] - weighted_mean[:, None]) ** 2, axis=1) / (k - 1)
    denominator = 1 + 2 * (k - 2) * correction / (k**2 - 1)
    statistic[valid] = numerator / denominator
    statistic = np.where(np.isfinite(statistic) & (statistic >= 0), statistic, 0.0)
    return statistic, df2


def _cluster_results(
    statistic: FloatArray,
    labels: IntArray,
    masses: FloatArray,
    null_maxima: FloatArray,
    alpha: float,
    time: FloatArray | None,
    *,
    f_statistic: bool,
) -> tuple[NDArray[np.bool_], FloatArray, pd.DataFrame]:
    """Return sample mask, cluster-level p map, and significant-only table."""
    significant = np.zeros(statistic.shape, dtype=bool)
    corrected = np.full(statistic.shape, np.nan)
    rows = []
    for cluster_id, mass in enumerate(masses, start=1):
        members = np.flatnonzero(labels == cluster_id)
        cluster_p = (1 + np.count_nonzero(null_maxima >= mass)) / (null_maxima.size + 1)
        corrected[members] = cluster_p
        if cluster_p > alpha:
            continue
        significant[members] = True
        cluster_stats = statistic[members]
        peak = float(np.max(cluster_stats)) if f_statistic else float(cluster_stats[np.argmax(np.abs(cluster_stats))])
        row = {
            "start": int(members[0]),
            "end": int(members[-1]),
            "MeanStatistic": float(np.mean(cluster_stats)),
            "PeakStatistic": peak,
            "Mass": float(mass),
            "CorrectedP": float(cluster_p),
        }
        rows.append(row)
    rows.sort(key=lambda row: row["start"])
    if time is None:
        columns = ["Cluster", "StartSample", "EndSample", "MeanStatistic", "PeakStatistic", "Mass", "CorrectedP"]
    else:
        columns = ["Cluster", "StartTime", "EndTime", "MeanStatistic", "PeakStatistic", "Mass", "CorrectedP"]
    table_rows = []
    for number, row in enumerate(rows, start=1):
        if time is None:
            bounds = {"StartSample": row["start"] + 1, "EndSample": row["end"] + 1}
        else:
            bounds = {"StartTime": float(time[row["start"]]), "EndTime": float(time[row["end"]])}
        table_rows.append({"Cluster": number, **bounds, **{key: row[key] for key in columns[3:]}})
    return significant, corrected, pd.DataFrame(table_rows, columns=columns)


__all__ = [
    "TemporalClusterSettings", "TemporalClusterStats", "_series", "_options",
    "_integer_vector", "_pairs", "_clusters", "_paired_t", "_repeated_f",
    "_welch_t", "_welch_f", "_cluster_results",
]
