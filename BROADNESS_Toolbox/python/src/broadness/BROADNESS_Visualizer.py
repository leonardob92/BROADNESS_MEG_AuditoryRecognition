"""Publication-ready visualization of BROADNESS results.

=========================================================================
 USER DOCUMENTATION
=========================================================================

This module provides the public function ``BROADNESS_Visualizer``. It takes
the output of ``BROADNESS_NetworkEstimation`` (or the future ICA equivalent)
and can produce the same five output families as the MATLAB visualizer:

1. Dynamic brain-activity maps of the original data.
2. Variance explained by the brain networks.
3. Network time series, including standard errors for participant data.
4. A three-dimensional view of spatial activation patterns in MNI space.
5. NIfTI images and Excel tables of thresholded spatial patterns.

OPTIONS
-------
Options are supplied in a dictionary using the established MATLAB names:

``WhichPlots``
    Five binary values selecting the outputs listed above. If omitted, all
    five outputs are requested.
``name_nii``
    Output directory used when NIfTI export is selected.
``MNI_coords``
    MNI coordinates with dimensions ``sources x 3``. Required for the 3D
    spatial plot and NIfTI export.
``ncomps``
    One-based component numbers used in all plots except the variance plot.
``ncomps_var``
    Number of components displayed in the variance plot. The default is the
    first 20, or all available components when fewer than 20 exist.
``Labels``
    One label for every experimental condition.
``color_PCs``
    Optional RGB colours for the selected components.
``color_conds``
    Optional RGB colours for the experimental conditions.

EXAMPLE
-------
Generate activity, variance, time-series, and spatial-pattern figures::

    from broadness import BROADNESS_Visualizer

    Options = {
        "WhichPlots": [1, 1, 1, 1, 0],
        "MNI_coords": MNI_coords,
        "ncomps": [1, 2, 3],
        "Labels": ["Condition 1", "Condition 2"],
    }
    VISUALIZATION = BROADNESS_Visualizer(BROADNESS, Options)

To export NIfTI images and Excel tables, select the fifth output and provide
``name_nii``. Files are written inside
``name_nii/BROADNESS_Output/BROADNESS_nifti``.

OUTPUT
------
The function returns a ``BROADNESSVisualization`` object containing the
generated Matplotlib figures and exported-file paths. Users do not construct
this object themselves. Use, for example::

    VISUALIZATION.variance_figure
    VISUALIZATION.time_series_figures
    VISUALIZATION.nifti_paths

INDEXING
--------
Component numbers supplied in ``Options["ncomps"]`` are one-based for
consistency with MATLAB. Conversion to Python's zero-based indexing occurs
internally.

If you use this toolbox, please cite the first BROADNESS paper:
Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
F., Testa, C., Vuust, P., Kringelbach, M. L., & Rosso, M. (2025).
BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in
Auditory Memory Networks. Advanced Science.
https://doi.org/10.1002/advs.202507878

Toolbox authors
---------------
Leonardo Bonetti - Center for Music in the Brain, Aarhus University; Centre
for Eudaimonia and Human Flourishing, Linacre College, University of Oxford -
leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
Mattia Rosso - Center for Music in the Brain, Aarhus University -
mattia.rosso@clin.au.dk
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TYPE_CHECKING
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from matplotlib.figure import Figure


# =========================================================================
#  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS)
#  VISUALIZER
# =========================================================================

FloatArray = NDArray[np.float64]
IntegerArray = NDArray[np.int64]


@dataclass(slots=True)
class BROADNESSVisualization:
    """Figures and files returned by :func:`BROADNESS_Visualizer`.

    Users normally inspect or save these outputs and do not construct this
    container directly.
    """

    dynamic_activity_figures: list[Figure]
    variance_figure: Figure | None
    time_series_figures: list[Figure]
    spatial_pattern_figure: Figure | None
    nifti_paths: list[Path]
    activation_table_paths: list[Path]


@dataclass(slots=True)
class _VisualizerInputs:
    data: FloatArray
    time: FloatArray
    activation_patterns: FloatArray
    time_series: FloatArray
    variance: FloatArray | None
    variance_permutations: FloatArray | None
    is_pca: bool
    has_participants: bool


# =========================================================================
#  INTERNAL HELPER FUNCTIONS
#  Users do not need to call the underscore-prefixed functions below.
# =========================================================================


def _field(
    container: Any,
    name: str,
    *,
    required: bool = True,
    default: Any = None,
) -> Any:
    """Read a field from a mapping or an object with MATLAB-like fields."""
    if isinstance(container, Mapping):
        if name in container:
            return container[name]
    elif hasattr(container, name):
        return getattr(container, name)

    if required:
        raise ValueError(f'"{name}" is required')
    return default


def _real_array(values: Any, name: str) -> FloatArray:
    """Return finite, real numeric values as a floating-point array."""
    array = np.asarray(values)
    if not np.issubdtype(array.dtype, np.number) or np.iscomplexobj(array):
        raise TypeError(f'"{name}" must contain real numeric values')
    array = np.asarray(array, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f'"{name}" cannot contain NaN or infinite values')
    return array


def _canonical_inputs(BROADNESS: Any) -> _VisualizerInputs:
    """Validate BROADNESS fields and add singleton condition dimensions."""
    data_original = _real_array(
        _field(BROADNESS, "OriginalData"),
        "BROADNESS.OriginalData",
    )
    if data_original.ndim not in (2, 3, 4):
        raise ValueError(
            '"BROADNESS.OriginalData" must have two, three, or four dimensions'
        )
    has_participants = data_original.ndim == 4 and data_original.shape[3] > 1
    if data_original.ndim == 2:
        data = data_original[:, :, np.newaxis, np.newaxis]
    elif data_original.ndim == 3:
        data = data_original[:, :, :, np.newaxis]
    else:
        data = data_original

    time = _real_array(_field(BROADNESS, "Time"), "BROADNESS.Time")
    if time.ndim > 2 or (time.ndim == 2 and 1 not in time.shape):
        raise ValueError('"BROADNESS.Time" must be a vector')
    time = time.reshape(-1)
    if time.size != data.shape[1]:
        raise ValueError(
            'The length of "BROADNESS.Time" must match the time dimension '
            'of "BROADNESS.OriginalData"'
        )

    activation_patterns = _real_array(
        _field(BROADNESS, "ActivationPatterns_BrainNetworks"),
        "BROADNESS.ActivationPatterns_BrainNetworks",
    )
    if activation_patterns.ndim != 2:
        raise ValueError(
            '"BROADNESS.ActivationPatterns_BrainNetworks" must have shape '
            "sources x components"
        )
    if activation_patterns.shape[0] != data.shape[0]:
        raise ValueError(
            'The source dimension of "ActivationPatterns_BrainNetworks" '
            'must match "OriginalData"'
        )

    time_series_original = _real_array(
        _field(BROADNESS, "TimeSeries_BrainNetworks"),
        "BROADNESS.TimeSeries_BrainNetworks",
    )
    if time_series_original.size == 0:
        time_series = np.empty(
            (time.size, 0, data.shape[2], data.shape[3]),
            dtype=float,
        )
    elif time_series_original.ndim == 2:
        time_series = time_series_original[:, :, np.newaxis, np.newaxis]
    elif time_series_original.ndim == 3:
        time_series = time_series_original[:, :, :, np.newaxis]
    elif time_series_original.ndim == 4:
        time_series = time_series_original
    else:
        raise ValueError(
            '"BROADNESS.TimeSeries_BrainNetworks" must have two, three, '
            "or four dimensions"
        )
    if time_series.shape[0] != time.size:
        raise ValueError(
            'The first dimension of "TimeSeries_BrainNetworks" must match '
            'the length of "Time"'
        )
    if time_series.shape[2:] != data.shape[2:]:
        raise ValueError(
            'The condition and participant dimensions of "TimeSeries_'
            'BrainNetworks" must match "OriginalData"'
        )

    variance_value = _field(
        BROADNESS,
        "Variance_BrainNetworks",
        required=False,
    )
    variance: FloatArray | None = None
    if variance_value is not None and np.asarray(variance_value).size > 0:
        variance = _real_array(
            variance_value,
            "BROADNESS.Variance_BrainNetworks",
        ).reshape(-1)
        if np.any(variance < 0):
            raise ValueError('"Variance_BrainNetworks" cannot be negative')

    variance_permutation_value = _field(
        BROADNESS,
        "VariancePermutations",
        required=False,
    )
    variance_permutations: FloatArray | None = None
    if (
        variance_permutation_value is not None
        and np.asarray(variance_permutation_value).size > 0
    ):
        variance_permutations = _real_array(
            variance_permutation_value,
            "BROADNESS.VariancePermutations",
        ).reshape(-1)

    return _VisualizerInputs(
        data=data,
        time=time,
        activation_patterns=activation_patterns,
        time_series=time_series,
        variance=variance,
        variance_permutations=variance_permutations,
        is_pca=variance is not None,
        has_participants=has_participants,
    )


def _which_plots(Options: Any) -> tuple[bool, bool, bool, bool, bool]:
    """Validate the five MATLAB-compatible plot-selection flags."""
    values = np.asarray(
        _field(
            Options,
            "WhichPlots",
            required=False,
            default=np.ones(5, dtype=int),
        )
    )
    if values.ndim != 1 or values.size != 5:
        raise ValueError('"Options.WhichPlots" must contain five binary values')
    if not np.all(np.isin(values, [0, 1, False, True])):
        raise ValueError('"Options.WhichPlots" must contain only zeros and ones')
    return tuple(bool(value) for value in values)  # type: ignore[return-value]


def _positive_integer(value: Any, name: str) -> int:
    """Validate a positive integer without accepting booleans."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value,
        (int, np.integer),
    ):
        raise TypeError(f'"{name}" must be a positive integer')
    if int(value) < 1:
        raise ValueError(f'"{name}" must be a positive integer')
    return int(value)


def _selected_components(
    BROADNESS: Any,
    Options: Any,
    inputs: _VisualizerInputs,
) -> IntegerArray:
    """Resolve one-based user selections and convert them to zero-based."""
    supplied = _field(Options, "ncomps", required=False)
    if supplied is None:
        significant = _field(
            BROADNESS,
            "Significant_BrainNetworks",
            required=False,
        )
        if inputs.is_pca and significant is not None and not isinstance(
            significant,
            str,
        ):
            one_based = np.asarray(significant)
        else:
            available = min(
                inputs.activation_patterns.shape[1],
                inputs.time_series.shape[1],
            )
            one_based = np.arange(1, min(5, available) + 1)
    else:
        one_based = np.asarray(supplied)

    if one_based.ndim != 1:
        raise ValueError('"Options.ncomps" must be a one-dimensional vector')
    if one_based.size == 0:
        return np.empty(0, dtype=np.int64)
    if not np.issubdtype(one_based.dtype, np.number) or np.iscomplexobj(one_based):
        raise TypeError('"Options.ncomps" must contain integer component numbers')
    numeric = np.asarray(one_based, dtype=float)
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.floor(numeric)):
        raise ValueError('"Options.ncomps" must contain integer component numbers')

    maximum = inputs.activation_patterns.shape[1]
    if inputs.time_series.shape[1] > 0:
        maximum = min(maximum, inputs.time_series.shape[1])
    if np.any(numeric < 1) or np.any(numeric > maximum):
        raise ValueError(
            '"Options.ncomps" contains numbers outside the available brain networks'
        )
    if np.unique(numeric).size != numeric.size:
        raise ValueError('"Options.ncomps" cannot contain duplicate components')
    return numeric.astype(np.int64) - 1


def _labels(Options: Any, number_conditions: int) -> list[str]:
    """Return supplied condition labels or concise defaults."""
    supplied = _field(Options, "Labels", required=False)
    if supplied is None:
        return [f"Condition {index}" for index in range(1, number_conditions + 1)]
    if isinstance(supplied, str):
        labels = [supplied]
    else:
        labels = [str(value) for value in supplied]
    if len(labels) != number_conditions:
        raise ValueError(
            '"Options.Labels" must contain one label for every condition'
        )
    return labels


def _matplotlib():
    """Import Matplotlib with a focused optional-dependency error."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "BROADNESS visualization requires matplotlib. "
            "Install broadness[visualize]."
        ) from error
    return plt


def _publication_style() -> dict[str, Any]:
    """Return a clean figure style using Helvetica Neue when available."""
    from matplotlib import font_manager

    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    font_family = "DejaVu Sans"
    for candidate in ("Helvetica Neue", "Helvetica", "Arial"):
        if candidate in available_fonts:
            font_family = candidate
            break

    return {
        "font.family": font_family,
        "font.size": 10.5,
        "axes.titlesize": 13,
        "axes.titleweight": 500,
        "axes.labelsize": 11,
        "axes.labelweight": 400,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9.5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }


def _default_colors(number: int, kind: str) -> FloatArray:
    """Create distinguishable component or condition colours."""
    plt = _matplotlib()
    if number == 0:
        return np.empty((0, 3), dtype=float)
    if kind == "components":
        positions = np.linspace(0.1, 0.9, number)
        return np.asarray(plt.get_cmap("cool")(positions)[:, :3]) * 0.78
    if number <= 10:
        return np.asarray(plt.get_cmap("tab10")(np.arange(number))[:, :3]) * 0.86
    positions = np.linspace(0.05, 0.95, number)
    return np.asarray(plt.get_cmap("turbo")(positions)[:, :3]) * 0.86


def _colors(
    Options: Any,
    field: str,
    number: int,
    kind: str,
) -> FloatArray:
    """Validate user colours and fall back when too few are provided."""
    supplied = _field(Options, field, required=False)
    if supplied is None:
        return _default_colors(number, kind)
    colors = _real_array(supplied, f"Options.{field}")
    if colors.ndim != 2 or colors.shape[1] not in (3, 4):
        raise ValueError(f'"Options.{field}" must have shape colours x 3 or 4')
    if np.any(colors < 0) or np.any(colors > 1):
        raise ValueError(f'"Options.{field}" values must be between zero and one')
    if colors.shape[0] < number:
        warnings.warn(
            f'Fewer colours than required were supplied in "Options.{field}"; '
            "using publication defaults instead",
            UserWarning,
            stacklevel=3,
        )
        return _default_colors(number, kind)
    return colors[:number, :3]


def _format_axis(axis: Any, *, grid_axis: str = "both") -> None:
    """Apply consistent publication-ready Cartesian-axis formatting."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.minorticks_on()
    axis.grid(True, which="major", axis=grid_axis, color="#D5D9DE", linewidth=0.7)
    axis.grid(False, which="minor", axis=grid_axis)
    axis.set_axisbelow(True)


def _dynamic_activity_plots(
    data: FloatArray,
    time: FloatArray,
    labels: list[str],
) -> list[Figure]:
    """Plot participant-averaged source activity for every condition."""
    plt = _matplotlib()
    group_data = np.mean(data, axis=3)
    maximum = float(np.max(np.abs(group_data)))
    if maximum == 0:
        maximum = 1.0
    figures: list[Figure] = []

    print("Generating dynamic brain activity plots")
    for condition, label in enumerate(labels):
        figure, axis = plt.subplots(
            figsize=(9.0, 5.2),
            constrained_layout=True,
        )
        image = axis.imshow(
            group_data[:, :, condition],
            aspect="auto",
            origin="upper",
            interpolation="nearest",
            extent=(
                float(time[0]),
                float(time[-1]),
                data.shape[0] + 0.5,
                0.5,
            ),
            cmap="RdBu_r",
            vmin=-maximum,
            vmax=maximum,
            rasterized=True,
        )
        colorbar = figure.colorbar(image, ax=axis, pad=0.02, fraction=0.045)
        colorbar.set_label("Amplitude")
        axis.set_title(f"Dynamic brain activity — {label}")
        axis.set_xlabel("Time (s)")
        axis.set_ylabel("Brain source")
        axis.grid(False)
        figures.append(figure)
    return figures


def _variance_plot(
    variance: FloatArray,
    variance_permutations: FloatArray | None,
    number_components: int,
) -> Figure:
    """Plot the empirical and optional randomized eigenspectra."""
    plt = _matplotlib()
    components = np.arange(1, number_components + 1)
    figure, axis = plt.subplots(
        figsize=(8.3, 5.1),
        constrained_layout=True,
    )
    axis.plot(
        components,
        variance[:number_components],
        color="#1F5A85",
        linewidth=2.1,
        marker="o",
        markersize=5.5,
        markerfacecolor="white",
        markeredgewidth=1.4,
        label="Data",
    )
    if variance_permutations is not None:
        if variance_permutations.size < number_components:
            raise ValueError(
                '"VariancePermutations" contains fewer values than requested '
                'by "Options.ncomps_var"'
            )
        axis.plot(
            components,
            variance_permutations[:number_components],
            color="#7A7F87",
            linewidth=1.8,
            linestyle="--",
            marker="s",
            markersize=4.7,
            markerfacecolor="white",
            markeredgewidth=1.2,
            label="Randomized data",
        )
    axis.set_xlim(0.7, number_components + 0.3)
    axis.set_xticks(components)
    axis.set_ylim(bottom=0)
    axis.set_title("Variance explained by principal components")
    axis.set_xlabel("Component")
    axis.set_ylabel("Variance explained (%)")
    _format_axis(axis)
    axis.legend(frameon=False, loc="upper right")
    return figure


def _time_series_plots(
    inputs: _VisualizerInputs,
    components: IntegerArray,
    labels: list[str],
    colors: FloatArray,
) -> list[Figure]:
    """Plot condition time courses and participant standard errors."""
    plt = _matplotlib()
    time_series_mean = np.mean(inputs.time_series, axis=3)
    time_series_sem: FloatArray | None = None
    if inputs.has_participants:
        time_series_sem = np.std(inputs.time_series, axis=3, ddof=1) / np.sqrt(
            inputs.time_series.shape[3]
        )

    figures: list[Figure] = []
    print("Generating time series plots for brain networks")
    for component in components:
        figure, axis = plt.subplots(
            figsize=(8.8, 5.0),
            constrained_layout=True,
        )
        for condition, label in enumerate(labels):
            mean = time_series_mean[:, component, condition]
            if time_series_sem is not None:
                sem = time_series_sem[:, component, condition]
                axis.fill_between(
                    inputs.time,
                    mean - sem,
                    mean + sem,
                    color=colors[condition],
                    alpha=0.18,
                    linewidth=0,
                )
            axis.plot(
                inputs.time,
                mean,
                color=colors[condition],
                linewidth=2.0,
                label=label,
            )

        title = f"Brain network {component + 1} — time series"
        if inputs.variance is not None and component < inputs.variance.size:
            title += f"  |  Variance = {inputs.variance[component]:.2f}%"
        axis.set_title(title)
        axis.set_xlabel("Time (s)")
        axis.set_ylabel("Component amplitude")
        axis.set_xlim(float(inputs.time[0]), float(inputs.time[-1]))
        _format_axis(axis)
        axis.legend(
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.13),
            ncol=min(5, len(labels)),
        )
        figures.append(figure)
    return figures


def _sample_standard_deviation(values: FloatArray) -> float:
    """Match MATLAB's default sample standard deviation."""
    if values.size < 2:
        return 0.0
    return float(np.std(values, ddof=1))


def _positive_threshold(values: FloatArray) -> FloatArray:
    """Select values at least one SD above the signed pattern mean."""
    threshold = float(np.mean(values)) + _sample_standard_deviation(values)
    return np.isfinite(values) & (values >= threshold)


def _absolute_threshold(values: FloatArray) -> FloatArray:
    """Threshold spatial magnitudes at mean absolute value plus one SD."""
    absolute = np.abs(values)
    threshold = float(np.mean(absolute)) + _sample_standard_deviation(absolute)
    return np.isfinite(values) & (absolute >= threshold)


def _equal_3d_axes(axis: Any, coordinates: FloatArray) -> None:
    """Apply equal physical scaling to MNI x, y, and z dimensions."""
    minimum = np.min(coordinates, axis=0)
    maximum = np.max(coordinates, axis=0)
    center = (minimum + maximum) / 2.0
    radius = max(float(np.max(maximum - minimum)) / 2.0, 1.0)
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1, 1, 1))


def _spatial_pattern_plot(
    activation_patterns: FloatArray,
    coordinates: FloatArray,
    components: IntegerArray,
    colors: FloatArray,
) -> Figure:
    """Plot selected thresholded patterns over a neutral MNI coordinate cloud."""
    plt = _matplotlib()
    from matplotlib.lines import Line2D

    valid_coordinates = np.all(np.isfinite(coordinates), axis=1)
    brain_coordinates = coordinates[valid_coordinates]
    if brain_coordinates.size == 0:
        raise ValueError('"Options.MNI_coords" contains no valid coordinates')

    figure = plt.figure(figsize=(9.2, 7.3), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    axis.scatter(
        brain_coordinates[:, 0],
        brain_coordinates[:, 1],
        brain_coordinates[:, 2],
        s=4,
        color="#59616B",
        alpha=0.075,
        linewidths=0,
        depthshade=False,
        rasterized=True,
    )

    legend_handles: list[Line2D] = []
    for position, component in enumerate(components):
        values = activation_patterns[:, component]
        mask = valid_coordinates & _positive_threshold(values)
        if np.any(mask):
            selected_values = values[mask]
            value_range = float(np.ptp(selected_values))
            if value_range > 0:
                normalized = (selected_values - np.min(selected_values)) / value_range
            else:
                normalized = np.ones(selected_values.size)
            order = np.argsort(normalized, kind="stable")
            selected_coordinates = coordinates[mask][order]
            marker_sizes = 14.0 + 82.0 * normalized[order]
            axis.scatter(
                selected_coordinates[:, 0],
                selected_coordinates[:, 1],
                selected_coordinates[:, 2],
                s=marker_sizes,
                color=colors[position],
                alpha=0.86,
                edgecolors="white",
                linewidths=0.35,
                depthshade=False,
                rasterized=True,
            )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=colors[position],
                markeredgecolor="white",
                markersize=9,
                label=f"Network {component + 1}",
            )
        )

    axis.view_init(elev=18, azim=-68)
    _equal_3d_axes(axis, brain_coordinates)
    axis.set_axis_off()
    axis.set_title("Spatial activation patterns of brain networks", pad=14)
    if legend_handles:
        axis.legend(
            handles=legend_handles,
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(0.0, 0.98),
        )
    return figure


def _mni_coordinates(Options: Any, number_sources: int) -> FloatArray:
    """Validate one MNI coordinate triplet per source."""
    coordinates_value = np.asarray(_field(Options, "MNI_coords"))
    if not np.issubdtype(coordinates_value.dtype, np.number) or np.iscomplexobj(
        coordinates_value
    ):
        raise TypeError('"Options.MNI_coords" must contain real numeric values')
    coordinates = np.asarray(coordinates_value, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape != (number_sources, 3):
        raise ValueError(
            '"Options.MNI_coords" must have shape sources x 3 and match '
            '"ActivationPatterns_BrainNetworks"'
        )
    if np.any(np.isinf(coordinates)):
        raise ValueError('"Options.MNI_coords" cannot contain infinite values')
    if not np.any(np.all(np.isfinite(coordinates), axis=1)):
        raise ValueError('"Options.MNI_coords" contains no valid coordinate rows')
    return coordinates


def _default_template_path() -> Path:
    """Locate the 8-mm template bundled with the MATLAB toolbox."""
    toolbox_root = Path(__file__).resolve().parents[3]
    return toolbox_root / "BROADNESS_External" / "MNI152_8mm_brain_diy.nii.gz"


def _nifti_and_table_output(
    activation_patterns: FloatArray,
    coordinates: FloatArray,
    components: IntegerArray,
    output_root: Path,
    *,
    is_pca: bool,
) -> tuple[list[Path], list[Path]]:
    """Export independent thresholded component volumes and activation tables."""
    try:
        import nibabel as nib
    except ImportError as error:
        raise ImportError(
            "NIfTI export requires nibabel. Install broadness[visualize]."
        ) from error
    try:
        import pandas as pd
    except ImportError as error:
        raise ImportError(
            "Activation-table export requires pandas and openpyxl. "
            "Install broadness[visualize]."
        ) from error

    template_path = _default_template_path()
    if not template_path.is_file():
        raise FileNotFoundError(f"NIfTI template not found: {template_path}")
    template = nib.load(str(template_path))
    shape = tuple(int(value) for value in template.shape[:3])
    inverse_affine = np.linalg.inv(template.affine)

    valid_coordinates = np.all(np.isfinite(coordinates), axis=1)
    voxel_indices = np.rint(
        nib.affines.apply_affine(inverse_affine, coordinates[valid_coordinates])
    ).astype(np.int64)
    in_bounds = np.all(voxel_indices >= 0, axis=1)
    in_bounds &= np.all(voxel_indices < np.asarray(shape), axis=1)
    coordinate_rows = np.flatnonzero(valid_coordinates)[in_bounds]
    voxel_indices = voxel_indices[in_bounds]

    output_directory = output_root / "BROADNESS_Output" / "BROADNESS_nifti"
    output_directory.mkdir(parents=True, exist_ok=True)
    nifti_paths: list[Path] = []
    table_paths: list[Path] = []
    method = "PCA" if is_pca else "ICA"

    print("Generating and saving NIfTI images of brain network activation patterns")
    for component in components:
        values = activation_patterns[:, component]
        threshold_mask = _absolute_threshold(values)
        thresholded = np.where(threshold_mask, values, 0.0)

        table_rows = np.flatnonzero(
            threshold_mask & valid_coordinates & (thresholded != 0)
        )
        if table_rows.size > 0:
            table = pd.DataFrame(
                {
                    "Index": np.arange(1, table_rows.size + 1),
                    "X": coordinates[table_rows, 0],
                    "Y": coordinates[table_rows, 1],
                    "Z": coordinates[table_rows, 2],
                    "Activation": thresholded[table_rows],
                }
            )
            table_path = output_directory / (
                f"ActivationTable_BrainNetwork_{component + 1}.xlsx"
            )
            try:
                table.to_excel(table_path, index=False)
            except ImportError as error:
                raise ImportError(
                    "Excel export requires openpyxl. "
                    "Install broadness[visualize]."
                ) from error
            print(f"Saved activation table: {table_path}")
            table_paths.append(table_path)
        else:
            print(
                f"No suprathreshold sources for brain network {component + 1}; "
                "skipping activation-table export"
            )

        # A fresh volume is essential so that components remain independent.
        volume = np.zeros(shape, dtype=np.float32)
        mapped_values = thresholded[coordinate_rows]
        nonzero = mapped_values != 0
        voxels = voxel_indices[nonzero]
        volume[voxels[:, 0], voxels[:, 1], voxels[:, 2]] = mapped_values[
            nonzero
        ].astype(np.float32)

        nifti_path = output_directory / (
            f"{method}_ActivationPattern_BrainNetwork_#{component + 1}.nii.gz"
        )
        image = nib.Nifti1Image(
            volume,
            template.affine,
            template.header.copy(),
        )
        image.set_data_dtype(np.float32)
        nib.save(image, str(nifti_path))
        print(f"Saved NIfTI image: {nifti_path}")
        nifti_paths.append(nifti_path)

    return nifti_paths, table_paths


# =========================================================================
#  MAIN PUBLIC FUNCTION
# =========================================================================


def BROADNESS_Visualizer(
    BROADNESS: Any,
    Options: Any,
    *,
    show: bool = True,
) -> BROADNESSVisualization:
    """Visualize the outputs of a BROADNESS PCA or ICA analysis.

    Parameters
    ----------
    BROADNESS
        Output from ``BROADNESS_NetworkEstimation`` or the corresponding ICA
        function. A mapping with the same MATLAB-style field names is also
        accepted.
    Options
        Mapping or object containing the MATLAB-compatible visualizer fields
        documented at the beginning of this module.
    show
        Call ``matplotlib.pyplot.show`` after preparing all figures. Set this
        to ``False`` for scripts, tests, or manual figure composition.

    Returns
    -------
    BROADNESSVisualization
        Generated figures and paths of exported NIfTI and Excel files.

    Notes
    -----
    Participant-level network time series are summarized using their mean and
    standard error. Spatial patterns remain group-level outputs of network
    estimation. Component selections in ``Options.ncomps`` are one-based.
    """
    if not isinstance(show, (bool, np.bool_)):
        raise TypeError('"show" must be a boolean')
    if Options is None:
        raise TypeError('"Options" must be a mapping or object with option fields')

    print("Checking visualizer inputs")
    inputs = _canonical_inputs(BROADNESS)
    which_plots = _which_plots(Options)
    components = _selected_components(BROADNESS, Options, inputs)
    labels = _labels(Options, inputs.data.shape[2])
    if (
        which_plots[2]
        and components.size > 0
        and inputs.time_series.shape[1] <= int(np.max(components))
    ):
        raise ValueError(
            "The selected network time series are unavailable, possibly because "
            "no components survived Monte Carlo simulations"
        )

    number_variance_components = min(
        20,
        inputs.variance.size if inputs.variance is not None else 0,
    )
    supplied_ncomps_var = _field(Options, "ncomps_var", required=False)
    if supplied_ncomps_var is not None:
        number_variance_components = _positive_integer(
            supplied_ncomps_var,
            "Options.ncomps_var",
        )
    if (
        inputs.variance is not None
        and number_variance_components > inputs.variance.size
    ):
        raise ValueError(
            '"Options.ncomps_var" exceeds the available principal components'
        )

    component_colors = _colors(
        Options,
        "color_PCs",
        components.size,
        "components",
    )
    condition_colors = _colors(
        Options,
        "color_conds",
        inputs.data.shape[2],
        "conditions",
    )

    coordinates: FloatArray | None = None
    if which_plots[3] or which_plots[4]:
        coordinates = _mni_coordinates(
            Options,
            inputs.activation_patterns.shape[0],
        )

    output_root: Path | None = None
    if which_plots[4]:
        name_nii = _field(Options, "name_nii")
        if not isinstance(name_nii, (str, Path)):
            raise TypeError('"Options.name_nii" must be a string or pathlib.Path')
        output_root = Path(name_nii).expanduser()

    plt = _matplotlib()
    dynamic_figures: list[Figure] = []
    variance_figure: Figure | None = None
    time_series_figures: list[Figure] = []
    spatial_pattern_figure: Figure | None = None
    nifti_paths: list[Path] = []
    table_paths: list[Path] = []

    with plt.rc_context(_publication_style()):
        if which_plots[0]:
            dynamic_figures = _dynamic_activity_plots(
                inputs.data,
                inputs.time,
                labels,
            )

        if which_plots[1]:
            if not inputs.is_pca:
                warnings.warn(
                    "Variance information is unavailable for this result; "
                    "skipping the variance plot",
                    UserWarning,
                    stacklevel=2,
                )
            else:
                print("Generating variance explained plot")
                variance_figure = _variance_plot(
                    inputs.variance,
                    inputs.variance_permutations,
                    number_variance_components,
                )

        if which_plots[2]:
            time_series_figures = _time_series_plots(
                inputs,
                components,
                labels,
                condition_colors,
            )

        if which_plots[3]:
            print("Generating 3D topographic plot of brain networks")
            spatial_pattern_figure = _spatial_pattern_plot(
                inputs.activation_patterns,
                coordinates,
                components,
                component_colors,
            )

        if which_plots[4]:
            nifti_paths, table_paths = _nifti_and_table_output(
                inputs.activation_patterns,
                coordinates,
                components,
                output_root,
                is_pca=inputs.is_pca,
            )

    if show:
        plt.show()

    return BROADNESSVisualization(
        dynamic_activity_figures=dynamic_figures,
        variance_figure=variance_figure,
        time_series_figures=time_series_figures,
        spatial_pattern_figure=spatial_pattern_figure,
        nifti_paths=nifti_paths,
        activation_table_paths=table_paths,
    )


__all__ = ["BROADNESSVisualization", "BROADNESS_Visualizer"]
