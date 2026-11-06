"""Phase-space trajectories, recurrence plots, and RQA measures."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy.spatial.distance import pdist, squareform

from .BROADNESS_FigureSettings import BROADNESS_FigureSettings


def _runs(bits: np.ndarray) -> np.ndarray:
    padded = np.r_[False, bits.astype(bool), False].astype(np.int8)
    starts = np.flatnonzero(np.diff(padded) == 1)
    ends = np.flatnonzero(np.diff(padded) == -1)
    return ends - starts


def _line_lengths(plot: np.ndarray, diagonal: bool) -> np.ndarray:
    if diagonal:
        lines = [np.diagonal(plot, offset=k) for k in range(1 - len(plot), len(plot))]
    else:
        lines = [plot[:, k] for k in range(plot.shape[1])]
    lengths = [_runs(line) for line in lines]
    return np.concatenate(lengths) if lengths else np.empty(0, dtype=int)


def _rqa_metrics(plot: np.ndarray, eligible: int, theiler_enabled: bool) -> list[float]:
    points = int(np.count_nonzero(plot))
    diagonal = np.sort(_line_lengths(plot, True))
    vertical = _line_lengths(plot, False)
    diagonal = diagonal[diagonal >= 2]
    vertical = vertical[vertical >= 2]
    entropy = 0.0
    if diagonal.size:
        frequencies = np.bincount(diagonal)[2:].astype(float)
        probabilities = frequencies[frequencies > 0] / diagonal.size
        entropy = float(-np.sum(probabilities * np.log(probabilities)))
    # MATLAB's historical DIV calculation drops the last (longest) sorted
    # diagonal when the main diagonal is retained.
    div_lines = diagonal if theiler_enabled or diagonal.size < 2 else diagonal[:-1]
    return [
        points / eligible,
        float(np.mean(diagonal)) if diagonal.size else np.nan,
        float(np.sum(diagonal) / points) if points else np.nan,
        entropy,
        float(np.mean(vertical)) if vertical.size else np.nan,
        float(np.sum(vertical) / points) if vertical.size and points else np.nan,
        float(np.max(vertical)) if vertical.size else np.nan,
        float(1 / np.max(div_lines)) if div_lines.size else np.nan,
    ]


def _field(value: object, name: str) -> object:
    if isinstance(value, dict):
        return value[name]
    return getattr(value, name)


def BROADNESS_PhaseSpace_RQA(
    BROADNESS: object,
    *,
    principalcomps: ArrayLike = (1, 2),
    timeinterval: ArrayLike | None = None,
    threshold: float = 0.1,
    normalization: str = "none",
    theiler_window: int | None = None,
    video: str = "off",
    figure: str = "off",
    figuremode: str | None = None,
    figurelayout: str = "individual",
    OutputPath: str | None = None,
    figureformats: str | tuple[str, ...] = ("png",),
    figureprefix: str = "",
    outpath: str | None = None,
) -> tuple[SimpleNamespace, SimpleNamespace]:
    """Analyse network scores (time × components × conditions × participants).

    Public component numbers are one-based. Recurrence plots use a strict
    distance threshold: distance < ``threshold * max(distance)``. The optional
    Theiler window removes pairs at or within that many samples of the main
    diagonal. Pairwise distance arrays require quadratic memory in time.
    """
    scores = np.asarray(_field(BROADNESS, "TimeSeries_BrainNetworks"), dtype=float)
    if scores.ndim < 2 or scores.ndim > 4 or not np.all(np.isfinite(scores)):
        raise ValueError("TimeSeries_BrainNetworks must be a finite 2D–4D array")
    scores = scores.reshape(scores.shape + (1,) * (4 - scores.ndim))
    times = np.asarray(_field(BROADNESS, "Time"), dtype=float).ravel()
    if times.size != scores.shape[0] or times.size < 2 or not np.all(np.diff(times) > 0):
        raise ValueError("Time must match the first dimension and increase strictly")
    pcs = np.asarray(principalcomps)
    if pcs.ndim != 1 or not pcs.size or not np.all(np.isfinite(pcs)) or np.any(pcs != np.floor(pcs)) or np.any((pcs < 1) | (pcs > scores.shape[1])):
        raise ValueError("principalcomps must contain available one-based components")
    pcs = pcs.astype(int)
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("threshold must be a finite nonnegative fraction")
    if normalization not in ("none", "pooled_zscore"):
        raise ValueError("normalization must be 'none' or 'pooled_zscore'")
    if theiler_window is not None and (not isinstance(theiler_window, (int, np.integer)) or theiler_window < 0):
        raise ValueError("theiler_window must be a nonnegative integer or None")
    if timeinterval is None:
        selected = np.arange(times.size)
    else:
        interval = np.asarray(timeinterval, dtype=float).ravel()
        if interval.size != 2 or not np.all(np.isfinite(interval)):
            raise ValueError("timeinterval must contain two finite times")
        lo, hi = sorted(interval)
        lo, hi = max(lo, times[0]), min(hi, times[-1])
        if lo >= hi:
            raise ValueError("timeinterval is outside the available time range")
        tolerance = max(10 * np.finfo(float).eps * max(1, np.max(np.abs(times))), np.median(np.diff(times)) * 1e-9)
        selected = np.flatnonzero((times >= lo - tolerance) & (times <= hi + tolerance))
    if selected.size < 2:
        raise ValueError("at least two time points are required")
    if theiler_window is not None and theiler_window >= selected.size - 1:
        raise ValueError("theiler_window must leave an eligible time pair")
    if figure not in ("on", "off") or video not in ("on", "off"):
        raise ValueError("figure and video must be 'on' or 'off'")
    mode = figuremode or ("show" if figure == "on" else "off")
    settings = BROADNESS_FigureSettings(mode, figurelayout,
                                        OutputPath if OutputPath is not None else outpath,
                                        figureformats, figureprefix, "PhaseSpace_RQA")

    trajectories = scores[np.ix_(selected, pcs - 1, np.arange(scores.shape[2]), np.arange(scores.shape[3]))].copy()
    center = np.zeros(pcs.size)
    scale = np.ones(pcs.size)
    if normalization == "pooled_zscore":
        center = trajectories.mean(axis=(0, 2, 3))
        scale = trajectories.std(axis=(0, 2, 3), ddof=1)
        if np.any(scale <= np.finfo(float).eps * np.maximum(1, np.max(np.abs(trajectories), axis=(0, 2, 3)))):
            raise ValueError("pooled z-score requires nonzero variance in each component")
        trajectories = (trajectories - center[None, :, None, None]) / scale[None, :, None, None]

    n_time, _, n_conditions, n_participants = trajectories.shape
    mask = np.ones((n_time, n_time), dtype=bool)
    if theiler_window is not None:
        indices = np.arange(n_time)
        mask = np.abs(indices[:, None] - indices[None, :]) > theiler_window
    distances = np.empty((n_conditions, n_participants), dtype=object)
    recurrences = np.empty_like(distances)
    tables: list[pd.DataFrame] = []
    columns = ["RR", "L", "DET", "ENTR", "TT", "LAM", "V_max", "DIV"]
    for participant in range(n_participants):
        rows = []
        for condition in range(n_conditions):
            coords = trajectories[:, :, condition, participant]
            distance = squareform(pdist(coords))
            plot = (distance < distance.max() * threshold) & mask
            distances[condition, participant] = distance
            recurrences[condition, participant] = plot.astype(int)
            rows.append(_rqa_metrics(plot, int(mask.sum()), theiler_window is not None))
        tables.append(pd.DataFrame(rows, columns=columns))

    participant_coordinates = np.empty((n_conditions, n_participants), dtype=object)
    mean_coordinates = np.empty((n_conditions, 1), dtype=object)
    for condition in range(n_conditions):
        mean_coordinates[condition, 0] = trajectories[:, :, condition, :].mean(axis=2)
        for participant in range(n_participants):
            participant_coordinates[condition, participant] = trajectories[:, :, condition, participant]
    figures = SimpleNamespace(Handles=[], Files=[], Animations=[])
    result = SimpleNamespace(
        PhaseSpace=SimpleNamespace(
            Time=times[selected], PCs=pcs, ParticipantCoordinates=participant_coordinates,
            MeanCoordinates=mean_coordinates, nConditions=n_conditions,
            nParticipants=n_participants, nDimensions=pcs.size,
            Normalization=SimpleNamespace(
                Method=normalization, Center=center, Scale=scale,
                Reference="Selected time-points pooled across all conditions and participants",
            ),
        ),
        RecurrencePlots=SimpleNamespace(DistMat=distances, RecurPlot=recurrences),
        RQA_metrics=tables,
        Figures=figures,
    )
    if settings.Mode != "off":
        from ._analysis_figures import rqa_figures

        figures = rqa_figures(result, threshold, theiler_window, settings)
        figures.Animations = []
        result.Figures = figures
    if video == "on":
        from ._analysis_figures import rqa_video

        animations, handles = rqa_video(result)
        figures.Animations.extend(animations)
        figures.Handles.extend(handles)
    return result, figures


__all__ = ["BROADNESS_PhaseSpace_RQA"]
