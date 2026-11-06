"""Time-resolved one-factor repeated-measures ANOVA with FDR correction."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import f as f_distribution

from .BROADNESS_FDRCorrection import BROADNESS_FDRCorrection, _probability


FloatArray = NDArray[np.float64]


@dataclass(slots=True)
class BROADNESSANOVAResult:
    """MATLAB-style fields from :func:`BROADNESS_RepeatedMeasuresANOVA`."""

    Enabled: bool
    FStatistics: FloatArray
    PValues: FloatArray
    UncorrectedPValues: FloatArray
    AdjustedPValues: FloatArray
    Significant: NDArray[np.bool_]
    FDRCriticalP: float
    Alpha: float
    FDRMethod: str
    FDRFamily: str
    GreenhouseGeisserEpsilon: FloatArray
    ConditionDF: FloatArray
    ErrorDF: FloatArray
    Correction: str


def BROADNESS_RepeatedMeasuresANOVA(
    data: ArrayLike, alpha: float = 0.05
) -> BROADNESSANOVAResult:
    """Test conditions at every time point and network.

    ``data`` has shape time × networks × conditions × participants. Complete
    participant rows are selected independently at each time/network cell.
    For more than two conditions, the p-values use Greenhouse–Geisser degrees
    of freedom. FDR is applied jointly across all tested cells.
    """
    level = _probability(alpha, "alpha")
    values = np.asarray(data)
    if values.ndim != 4 or not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise ValueError("data must be a real time x networks x conditions x participants array")
    values = np.asarray(values, dtype=float)
    n_time, n_networks, n_conditions, _ = values.shape
    if n_conditions < 2:
        raise ValueError("repeated-measures ANOVA requires at least two conditions")

    shape = (n_time, n_networks)
    f_statistics = np.full(shape, np.nan)
    p_values = np.full(shape, np.nan)
    uncorrected = np.full(shape, np.nan)
    epsilon_values = np.full(shape, np.nan)
    condition_df = np.full(shape, np.nan)
    error_df = np.full(shape, np.nan)

    for sample in range(n_time):
        for network in range(n_networks):
            cell = values[sample, network].T
            cell = cell[np.all(np.isfinite(cell), axis=1)]
            n_participants = cell.shape[0]
            if n_participants < 2:
                continue

            grand_mean = float(np.mean(cell))
            condition_means = np.mean(cell, axis=0)
            participant_means = np.mean(cell, axis=1)
            ss_condition = n_participants * np.sum((condition_means - grand_mean) ** 2)
            ss_participant = n_conditions * np.sum((participant_means - grand_mean) ** 2)
            ss_total = np.sum((cell - grand_mean) ** 2)
            ss_error = max(0.0, ss_total - ss_condition - ss_participant)
            df_condition = n_conditions - 1
            df_error = (n_participants - 1) * df_condition
            ms_condition = ss_condition / df_condition
            ms_error = ss_error / df_error
            if ms_error == 0:
                f_value = float("inf") if ms_condition > 0 else float("nan")
            else:
                f_value = float(ms_condition / ms_error)

            epsilon = 1.0
            if n_conditions > 2:
                centering = np.eye(n_conditions) - np.ones((n_conditions, n_conditions)) / n_conditions
                covariance = np.cov(cell, rowvar=False)
                centered_covariance = centering @ covariance @ centering
                denominator = df_condition * np.sum(centered_covariance ** 2)
                if denominator > 0:
                    epsilon = float(np.trace(centered_covariance) ** 2 / denominator)
                    epsilon = min(1.0, max(1.0 / df_condition, epsilon))

            f_statistics[sample, network] = f_value
            epsilon_values[sample, network] = epsilon
            condition_df[sample, network] = epsilon * df_condition
            error_df[sample, network] = epsilon * df_error
            uncorrected[sample, network] = f_distribution.sf(f_value, df_condition, df_error)
            p_values[sample, network] = f_distribution.sf(
                f_value, epsilon * df_condition, epsilon * df_error
            )

    significant, adjusted, critical = BROADNESS_FDRCorrection(p_values, level)
    return BROADNESSANOVAResult(
        Enabled=True,
        FStatistics=f_statistics,
        PValues=p_values,
        UncorrectedPValues=uncorrected,
        AdjustedPValues=adjusted,
        Significant=significant,
        FDRCriticalP=critical,
        Alpha=level,
        FDRMethod="Benjamini-Hochberg",
        FDRFamily="All tested time-points and networks",
        GreenhouseGeisserEpsilon=epsilon_values,
        ConditionDF=condition_df,
        ErrorDF=error_df,
        Correction="Greenhouse-Geisser for more than two conditions",
    )


__all__ = ["BROADNESSANOVAResult", "BROADNESS_RepeatedMeasuresANOVA"]
