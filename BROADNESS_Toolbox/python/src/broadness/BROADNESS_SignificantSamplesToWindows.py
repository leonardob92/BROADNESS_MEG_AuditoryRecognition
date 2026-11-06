"""Convert runs of significant time samples to plotting windows."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def BROADNESS_SignificantSamplesToWindows(
    significantSamples: ArrayLike,
    time: ArrayLike,
    minimumSamples: int = 1,
) -> list[tuple[float, float]]:
    """Return ``(start_time, end_time)`` for each retained contiguous run."""
    samples = np.asarray(significantSamples)
    times = np.asarray(time)
    samples_vector = samples.ndim == 1 or (samples.ndim == 2 and 1 in samples.shape)
    time_vector = times.ndim == 1 or (times.ndim == 2 and 1 in times.shape)
    if not samples_vector or not time_vector or samples.size != times.size:
        raise ValueError("significantSamples and time must be vectors of equal length")
    if not np.issubdtype(times.dtype, np.number) or np.iscomplexobj(times):
        raise TypeError("time must be a real numeric vector")
    samples = samples.reshape(-1)
    times = np.asarray(times, dtype=float).reshape(-1)
    if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError("time must contain finite, strictly increasing values")
    if isinstance(minimumSamples, (bool, np.bool_)) or not isinstance(
        minimumSamples, (int, np.integer)
    ) or minimumSamples < 1:
        raise ValueError("minimumSamples must be a positive integer")
    if not (np.issubdtype(samples.dtype, np.number) or samples.dtype == np.bool_):
        raise TypeError("significantSamples must be numeric or boolean")
    active = np.asarray(samples, dtype=bool)
    changes = np.diff(np.concatenate(([False], active, [False])).astype(int))
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1) - 1
    return [
        (float(times[start]), float(times[end]))
        for start, end in zip(starts, ends)
        if end - start + 1 >= minimumSamples
    ]


__all__ = ["BROADNESS_SignificantSamplesToWindows"]
