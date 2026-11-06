"""Benjamini–Hochberg correction over a user-defined family of p-values."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _probability(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise TypeError(f"{name} must be a real number between 0 and 1")
    result = float(value)
    if not np.isfinite(result) or not 0 < result < 1:
        raise ValueError(f"{name} must be strictly between 0 and 1")
    return result


def BROADNESS_FDRCorrection(
    pValues: ArrayLike, alpha: float = 0.05
) -> tuple[NDArray[np.bool_], NDArray[np.float64], float]:
    """Correct all finite p-values together, retaining NaNs in their positions.

    Returns ``(significant, adjusted_p, critical_p)``. ``critical_p`` is NaN
    when no hypothesis passes. The caller defines the correction family by
    supplying its p-values in one array.
    """
    level = _probability(alpha, "alpha")
    values = np.asarray(pValues)
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("pValues must be a real numeric array")
    values = np.asarray(values, dtype=float)
    if np.any((~np.isnan(values)) & ((values < 0) | (values > 1))):
        raise ValueError("finite p-values must be between 0 and 1")

    flat = values.ravel()
    valid = ~np.isnan(flat)
    adjusted = np.full(flat.shape, np.nan)
    significant = np.zeros(flat.shape, dtype=bool)
    critical = float("nan")
    if not np.any(valid):
        return significant.reshape(values.shape), adjusted.reshape(values.shape), critical

    valid_values = flat[valid]
    order = np.argsort(valid_values, kind="stable")
    sorted_p = valid_values[order]
    n_tests = sorted_p.size
    ranks = np.arange(1, n_tests + 1)
    passing = np.flatnonzero(sorted_p <= level * ranks / n_tests)
    if passing.size:
        critical = float(sorted_p[passing[-1]])
        significant[valid] = valid_values <= critical

    adjusted_sorted = np.minimum.accumulate(
        (sorted_p * n_tests / ranks)[::-1]
    )[::-1]
    adjusted_valid = np.empty(n_tests, dtype=float)
    adjusted_valid[order] = np.minimum(adjusted_sorted, 1.0)
    adjusted[valid] = adjusted_valid
    return significant.reshape(values.shape), adjusted.reshape(values.shape), critical


__all__ = ["BROADNESS_FDRCorrection"]
