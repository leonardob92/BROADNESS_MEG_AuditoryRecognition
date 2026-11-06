"""Participant-wise behavioural correlations with family-wide FDR correction."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import pearsonr, spearmanr

from .BROADNESS_FDRCorrection import BROADNESS_FDRCorrection


def BROADNESS_BehaviourCorrelation(
    data: ArrayLike,
    behaviour: ArrayLike,
    alpha: float = 0.05,
    correlationType: str = "Pearson",
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_], float]:
    """Correlate each measure with behaviour across the final participant axis.

    Missing observations are excluded separately for each measure. Results
    have the shape of ``data`` without its last axis. All resulting p-values
    are corrected together using Benjamini–Hochberg FDR.
    """
    values = np.asarray(data)
    if values.ndim < 2 or not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise ValueError("data must be a real numeric array with participants last")
    observations = np.asarray(values, dtype=float)
    behaviour_values = np.asarray(behaviour)
    if (
        behaviour_values.ndim not in (1, 2)
        or (behaviour_values.ndim == 2 and 1 not in behaviour_values.shape)
        or not np.issubdtype(behaviour_values.dtype, np.number)
        or np.iscomplexobj(behaviour_values)
        or behaviour_values.size != observations.shape[-1]
    ):
        raise ValueError("behaviour must have one real numeric value per participant")
    behaviour_values = np.asarray(behaviour_values, dtype=float).reshape(-1)
    if np.count_nonzero(np.isfinite(behaviour_values)) < 3:
        raise ValueError("behaviour needs at least three finite participant values")
    if not isinstance(correlationType, str) or correlationType.lower() not in (
        "pearson", "spearman"
    ):
        raise ValueError("correlationType must be Pearson or Spearman")

    output_shape = observations.shape[:-1]
    rows = observations.reshape(-1, observations.shape[-1])
    correlations = np.full(rows.shape[0], np.nan)
    p_values = np.full(rows.shape[0], np.nan)
    method = pearsonr if correlationType.lower() == "pearson" else spearmanr
    for index, row in enumerate(rows):
        valid = np.isfinite(row) & np.isfinite(behaviour_values)
        if np.count_nonzero(valid) < 2:
            continue
        if np.unique(row[valid]).size < 2 or np.unique(behaviour_values[valid]).size < 2:
            continue
        result = method(row[valid], behaviour_values[valid])
        correlations[index] = float(result.statistic)
        p_values[index] = float(result.pvalue)

    correlations = correlations.reshape(output_shape)
    p_values = p_values.reshape(output_shape)
    significant, adjusted, critical = BROADNESS_FDRCorrection(p_values, alpha)
    return correlations, p_values, adjusted, significant, critical


__all__ = ["BROADNESS_BehaviourCorrelation"]
