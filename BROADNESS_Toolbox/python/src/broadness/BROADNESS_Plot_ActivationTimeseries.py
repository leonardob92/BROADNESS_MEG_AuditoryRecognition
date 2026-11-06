"""Plot group and condition network time series with standard errors."""

from __future__ import annotations

from math import ceil, sqrt
from types import SimpleNamespace as Result
import warnings

import numpy as np
from numpy.typing import ArrayLike

from .BROADNESS_FigureSettings import BROADNESS_FigureSettings
from .BROADNESS_FinalizeFigure import BROADNESS_FinalizeFigure


def _limits(value: ArrayLike | None, name: str) -> tuple[float, float] | None:
    if value is None:
        return None
    limits = np.asarray(value, dtype=float).ravel()
    if limits.size != 2 or not np.all(np.isfinite(limits)) or limits[0] >= limits[1]:
        raise ValueError(f"{name} must be two increasing finite values")
    return float(limits[0]), float(limits[1])


def _windows(value: object, time: np.ndarray) -> list[list[tuple[float, float]]]:
    if value is None:
        return []
    result = []
    for level in value:
        windows = []
        for item in level:
            interval = np.asarray(item, dtype=float).ravel()
            if interval.size != 2 or not np.all(np.isfinite(interval)) or interval[0] > interval[1]:
                raise ValueError("each significant window must be [start, end]")
            if interval[0] < time[0] or interval[1] > time[-1]:
                warnings.warn("a significance window is outside the time range", UserWarning, stacklevel=3)
            windows.append((float(interval[0]), float(interval[1])))
        result.append(windows)
    return result


def _plot_network(
    axis: object, network: int, data: np.ndarray, time: np.ndarray,
    groups: list[np.ndarray], group_labels: list[str], condition_labels: list[str],
    colors: np.ndarray, ste_style: int, transparency: float,
    windows: list[list[tuple[float, float]]], significance_colors: np.ndarray | None,
    significance_labels: list[str] | None, significance_style: str,
    significance_position: str, significance_offset: float,
    significance_padding: float, xlimits: tuple[float, float] | None,
    ylimits: tuple[float, float] | None,
) -> None:
    all_values = []
    for condition in range(data.shape[2]):
        for group, indices in enumerate(groups):
            values = data[:, network, condition, :][:, indices]
            valid = np.isfinite(values)
            count = valid.sum(axis=1)
            mean = np.divide(np.where(valid, values, 0).sum(axis=1), count,
                             out=np.full(len(time), np.nan), where=count > 0)
            centered = np.where(valid, values - mean[:, None], 0)
            sd = np.sqrt(np.divide(np.sum(centered**2, axis=1), count - 1,
                                   out=np.zeros(len(time)), where=count > 1))
            se = sd / np.sqrt(len(indices))
            color = colors[condition, group]
            label = condition_labels[condition] if len(groups) == 1 else f"{condition_labels[condition]} - {group_labels[group]}"
            axis.plot(time, mean, color=color, linewidth=1.5, label=label, zorder=3)
            if ste_style == 1:
                axis.plot(time, mean + se, ":", color=color, linewidth=0.8, zorder=2)
                axis.plot(time, mean - se, ":", color=color, linewidth=0.8, zorder=2)
            else:
                axis.fill_between(time, mean - se, mean + se, color=color,
                                  alpha=transparency, linewidth=0, zorder=2)
            all_values.extend((mean + se, mean - se))
    finite = np.concatenate(all_values)
    finite = finite[np.isfinite(finite)]
    if ylimits is None:
        low, high = (float(finite.min()), float(finite.max())) if finite.size else (-1.0, 1.0)
        margin = max((high - low) * 0.05, np.finfo(float).eps)
        ylimits = (low - margin, high + margin)
    axis.set_ylim(*ylimits)
    if windows:
        n_levels = len(windows)
        span = ylimits[1] - ylimits[0]
        offset = significance_offset * span
        padding = significance_padding * span
        total = (n_levels - 1) * offset + padding
        if significance_position == "above":
            base = ylimits[1] + offset
            axis.set_ylim(ylimits[0], base + total)
            positions = [base + i * offset for i in range(n_levels)]
        else:
            base = ylimits[0] - offset
            axis.set_ylim(base - total, ylimits[1])
            positions = [base - i * offset for i in range(n_levels)]
        for level, intervals in enumerate(windows):
            default_color = (0.45, 0.45, 0.45) if significance_style == "line" else (0.85, 0.85, 0.85)
            color = default_color if significance_colors is None or level >= len(significance_colors) else significance_colors[level]
            for start, end in intervals:
                if significance_style == "patch":
                    axis.axvspan(start, end, color=color, alpha=0.2, zorder=0)
                else:
                    axis.plot((start, end), (positions[level], positions[level]),
                              color=color, linewidth=2, zorder=4)
            if significance_style == "line" and significance_labels and intervals:
                x_max = (xlimits or (float(time[0]), float(time[-1])))[1]
                x_min = (xlimits or (float(time[0]), float(time[-1])))[0]
                axis.text(x_max - 0.005 * (x_max - x_min), positions[level],
                          significance_labels[level], ha="right", va="center",
                          color=color, fontsize=8, fontweight="bold",
                          bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5})
    axis.set_xlim(*(xlimits or (float(time[0]), float(time[-1]))))
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Amplitude")
    axis.grid(True, alpha=0.2)
    axis.legend(loc="best", fontsize="small")


def BROADNESS_Plot_ActivationTimeseries(
    data: ArrayLike, time: ArrayLike, *, Groups: list[ArrayLike] | None = None,
    GroupLabels: list[str] | None = None, ConditionLabels: list[str] | None = None,
    ConditionXGroupColors: list[ArrayLike] | None = None,
    SignificantWindows: list[list[ArrayLike]] | None = None,
    SignificanceColors: ArrayLike | None = None,
    SignificanceLabels: list[str] | None = None,
    SignificanceStyle: str = "patch", SignificanceLinePosition: str = "above",
    SignificanceLineOffset: float = 0.01, SignificanceLinePadding: float = 0.05,
    YLimits: ArrayLike | None = None, XLimits: ArrayLike | None = None,
    STEStyle: int = 2, Transparency: float = 0.3,
    FigureMode: str = "show", FigureLayout: str = "individual",
    OutputPath: str | None = None, FigureFormats: str | tuple[str, ...] = ("png",),
    FigurePrefix: str = "",
) -> Result:
    """Plot data in time × networks × conditions × participants order.

    The older networks × time × participants × conditions order is accepted
    when the second dimension matches ``time``. Group member numbers are
    one-based for consistency with MATLAB.
    """
    settings = BROADNESS_FigureSettings(FigureMode, FigureLayout, OutputPath,
                                        FigureFormats, FigurePrefix, "ActivationTimeseries")
    values = np.asarray(data, dtype=float)
    times = np.asarray(time, dtype=float).ravel()
    if values.ndim < 2 or values.ndim > 4 or times.size < 2 or not np.all(np.diff(times) > 0):
        raise ValueError("data must be a 2D–4D array and time must increase")
    if values.shape[0] == len(times):
        values = values.reshape(values.shape + (1,) * (4 - values.ndim))
    elif values.shape[1] == len(times) and values.ndim == 4:
        values = values.transpose(1, 0, 3, 2)
    else:
        raise ValueError("time length must match the first or legacy second data dimension")
    n_time, n_networks, n_conditions, n_participants = values.shape
    if np.isnan(values).mean() > 0.05:
        warnings.warn("more than 5% of time-series data are NaN", UserWarning, stacklevel=2)
    groups = [np.arange(n_participants)]
    if Groups is not None:
        groups = []
        for row in Groups:
            indices = np.asarray(row)
            if not np.issubdtype(indices.dtype, np.number) or np.iscomplexobj(indices):
                raise ValueError("Groups must contain one-based integer participant indices")
            indices = np.asarray(indices, dtype=float).ravel()
            if np.any(~np.isfinite(indices)) or np.any(indices != np.floor(indices)):
                raise ValueError("Groups must contain one-based integer participant indices")
            groups.append(indices.astype(int) - 1)
    if not groups or any(len(row) == 0 or np.any((row < 0) | (row >= n_participants)) for row in groups):
        raise ValueError("Groups must contain nonempty one-based participant lists")
    group_labels = [f"Group {i}" for i in range(1, len(groups) + 1)] if GroupLabels is None else [str(x) for x in GroupLabels]
    condition_labels = [f"Condition {i}" for i in range(1, n_conditions + 1)] if ConditionLabels is None else [str(x) for x in ConditionLabels]
    if len(group_labels) != len(groups) or len(condition_labels) != n_conditions:
        raise ValueError("group and condition labels must match their counts")
    if ConditionXGroupColors is None:
        from matplotlib import pyplot as plt

        base = plt.get_cmap("tab10")
        colors = np.asarray([[base(c % 10)[:3] for _ in groups] for c in range(n_conditions)])
        for group in range(1, len(groups)):
            colors[:, group] *= max(0.3, 1 - 0.2 * group)
    else:
        colors = np.asarray(ConditionXGroupColors, dtype=float)
        if colors.shape != (n_conditions, len(groups), 3) or np.any((colors < 0) | (colors > 1)):
            raise ValueError("ConditionXGroupColors must be conditions × groups × RGB")
    if STEStyle not in (1, 2) or not 0 <= Transparency <= 1:
        raise ValueError("STEStyle must be 1 or 2 and Transparency must be in [0,1]")
    if SignificanceStyle not in ("patch", "line") or SignificanceLinePosition not in ("above", "below"):
        raise ValueError("invalid significance style or line position")
    if SignificanceLineOffset < 0 or SignificanceLinePadding < 0:
        raise ValueError("significance line offset and padding must be nonnegative")
    windows = _windows(SignificantWindows, times)
    significance_colors = None if SignificanceColors is None else np.asarray(SignificanceColors, dtype=float)
    if significance_colors is not None and (significance_colors.ndim != 2 or significance_colors.shape[1] != 3 or np.any((significance_colors < 0) | (significance_colors > 1))):
        raise ValueError("SignificanceColors must be levels × RGB")
    significance_labels = None if SignificanceLabels is None else [str(x) for x in SignificanceLabels]
    if significance_labels is not None and len(significance_labels) != len(windows):
        raise ValueError("SignificanceLabels must match SignificantWindows levels")
    xlimits, ylimits = _limits(XLimits, "XLimits"), _limits(YLimits, "YLimits")
    figures = Result(Handles=[], Files=[])
    if settings.Mode == "off":
        return figures
    from matplotlib import pyplot as plt

    args = (values, times, groups, group_labels, condition_labels, colors,
            STEStyle, Transparency, windows, significance_colors,
            significance_labels, SignificanceStyle, SignificanceLinePosition,
            SignificanceLineOffset, SignificanceLinePadding, xlimits, ylimits)
    if settings.MakeIndividual:
        for network in range(n_networks):
            fig, axis = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
            _plot_network(axis, network, *args)
            axis.set_title(f"Network {network + 1}")
            figures.Files.extend(BROADNESS_FinalizeFigure(fig, settings, f"Network_{network + 1:02d}"))
            if settings.Show:
                figures.Handles.append(fig)
    if settings.MakeSummary:
        columns = min(3, ceil(sqrt(n_networks)))
        rows = ceil(n_networks / columns)
        fig, axes = plt.subplots(rows, columns, figsize=(5 * columns, 3.8 * rows),
                                 constrained_layout=True, squeeze=False)
        fig.suptitle("Brain Network Activation Time Series")
        for network, axis in enumerate(axes.flat):
            if network < n_networks:
                _plot_network(axis, network, *args)
                axis.set_title(f"Network {network + 1}")
            else:
                axis.set_visible(False)
        figures.Files.extend(BROADNESS_FinalizeFigure(fig, settings, "Networks_Summary"))
        if settings.Show:
            figures.Handles.append(fig)
    return figures


__all__ = ["BROADNESS_Plot_ActivationTimeseries"]
