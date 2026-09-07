from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.io import loadmat

from broadness import (
    BROADNESSResult,
    BROADNESSVisualization,
    BROADNESS_Visualizer,
)


TOOLBOX_ROOT = Path(__file__).resolve().parents[2]
COORDINATE_FILE = (
    TOOLBOX_ROOT / "BROADNESS_External" / "MNI152_8mm_coord_dyi.mat"
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _visualizer_inputs() -> tuple[BROADNESSResult, np.ndarray]:
    generator = np.random.default_rng(21)
    sources = 8
    time = np.linspace(-0.2, 0.8, 51)
    conditions = 3
    participants = 4

    patterns, _ = np.linalg.qr(generator.normal(size=(sources, sources)))
    source = np.arange(1, sources + 1)[:, np.newaxis, np.newaxis, np.newaxis]
    sample = np.arange(time.size)[np.newaxis, :, np.newaxis, np.newaxis]
    condition = np.arange(1, conditions + 1)[np.newaxis, np.newaxis, :, np.newaxis]
    participant = np.arange(1, participants + 1)[np.newaxis, np.newaxis, np.newaxis, :]
    data = (
        np.sin(0.07 * source * sample)
        + 0.15 * np.cos(0.05 * (source + condition) * sample)
        + 0.02 * participant * np.sin(0.11 * sample)
    )
    time_series = np.einsum("stcp,sk->tkcp", data, patterns)
    variance = np.array([52, 20, 11, 7, 4, 3, 2, 1], dtype=float)
    coordinates = loadmat(COORDINATE_FILE)["MNI8"][:sources].astype(float)

    result = BROADNESSResult(
        Variance_BrainNetworks=variance,
        Significant_BrainNetworks="MCS has not been run",
        ActivationPatterns_BrainNetworks=patterns,
        TimeSeries_BrainNetworks=time_series,
        Time=time,
        OriginalData=data,
        VariancePermutations=None,
    )
    return result, coordinates


def test_publication_figures_cover_the_first_four_plot_families():
    result, coordinates = _visualizer_inputs()
    options = {
        "WhichPlots": [1, 1, 1, 1, 0],
        "MNI_coords": coordinates,
        "ncomps": [1, 3],
        "ncomps_var": 5,
        "Labels": ["Old", "New A", "New B"],
    }

    visualization = BROADNESS_Visualizer(result, options, show=False)

    assert isinstance(visualization, BROADNESSVisualization)
    assert len(visualization.dynamic_activity_figures) == 3
    assert visualization.variance_figure is not None
    assert len(visualization.time_series_figures) == 2
    assert visualization.spatial_pattern_figure is not None
    spatial_axis = visualization.spatial_pattern_figure.axes[0]
    assert any(
        isinstance(collection, Poly3DCollection)
        for collection in spatial_axis.collections
    )
    assert visualization.nifti_paths == []
    assert visualization.activation_table_paths == []

    time_axis = visualization.time_series_figures[0].axes[0]
    assert "Brain network 1" in time_axis.get_title()
    assert len(time_axis.lines) == 3
    assert len(time_axis.collections) == 3
    assert time_axis.get_xlabel() == "Time (s)"


def test_dynamic_activity_plots_share_a_symmetric_color_scale():
    result, _ = _visualizer_inputs()

    visualization = BROADNESS_Visualizer(
        result,
        {"WhichPlots": [1, 0, 0, 0, 0]},
        show=False,
    )

    limits = [
        figure.axes[0].images[0].get_clim()
        for figure in visualization.dynamic_activity_figures
    ]
    assert len(set(limits)) == 1
    assert limits[0][0] == -limits[0][1]


@pytest.mark.parametrize("number_dimensions", [2, 3])
def test_two_and_three_dimensional_results_are_supported(number_dimensions):
    result, _ = _visualizer_inputs()
    if number_dimensions == 2:
        result.OriginalData = result.OriginalData[:, :, 0, 0]
        result.TimeSeries_BrainNetworks = result.TimeSeries_BrainNetworks[:, :, 0, 0]
    else:
        result.OriginalData = result.OriginalData[:, :, :, 0]
        result.TimeSeries_BrainNetworks = result.TimeSeries_BrainNetworks[:, :, :, 0]

    visualization = BROADNESS_Visualizer(
        result,
        {"WhichPlots": [1, 0, 1, 0, 0], "ncomps": [1]},
        show=False,
    )

    expected_conditions = 1 if number_dimensions == 2 else 3
    assert len(visualization.dynamic_activity_figures) == expected_conditions
    assert len(visualization.time_series_figures) == 1


def test_variance_plot_includes_randomized_spectrum():
    result, _ = _visualizer_inputs()
    result.VariancePermutations = np.linspace(8, 1, 8)
    result.Significant_BrainNetworks = np.array([1, 2])

    visualization = BROADNESS_Visualizer(
        result,
        {"WhichPlots": [0, 1, 0, 0, 0], "ncomps_var": 6},
        show=False,
    )

    assert visualization.variance_figure is not None
    assert len(visualization.variance_figure.axes[0].lines) == 2


def test_ica_like_results_skip_only_the_variance_plot():
    result, _ = _visualizer_inputs()
    ica_result = SimpleNamespace(
        OriginalData=result.OriginalData,
        Time=result.Time,
        ActivationPatterns_BrainNetworks=result.ActivationPatterns_BrainNetworks,
        TimeSeries_BrainNetworks=result.TimeSeries_BrainNetworks,
    )

    with pytest.warns(UserWarning, match="Variance information"):
        visualization = BROADNESS_Visualizer(
            ica_result,
            {
                "WhichPlots": [0, 1, 1, 0, 0],
                "ncomps": [2],
                "Labels": ["A", "B", "C"],
            },
            show=False,
        )

    assert visualization.variance_figure is None
    assert len(visualization.time_series_figures) == 1
    assert "Variance" not in visualization.time_series_figures[0].axes[0].get_title()


def test_nifti_and_excel_outputs_are_independent_between_components(tmp_path):
    result, coordinates = _visualizer_inputs()
    patterns = np.zeros_like(result.ActivationPatterns_BrainNetworks)
    patterns[0, 0] = 10
    patterns[1, 1] = -12
    result.ActivationPatterns_BrainNetworks = patterns

    visualization = BROADNESS_Visualizer(
        result,
        {
            "WhichPlots": [0, 0, 0, 0, 1],
            "MNI_coords": coordinates,
            "name_nii": tmp_path,
            "ncomps": [1, 2],
        },
        show=False,
    )

    assert len(visualization.nifti_paths) == 2
    assert len(visualization.activation_table_paths) == 2
    assert all(path.is_file() for path in visualization.nifti_paths)
    assert all(path.is_file() for path in visualization.activation_table_paths)

    volumes = [np.asarray(nib.load(path).dataobj) for path in visualization.nifti_paths]
    assert np.count_nonzero(volumes[0]) == 1
    assert np.count_nonzero(volumes[1]) == 1
    assert not np.array_equal(volumes[0], volumes[1])

    first_table = pd.read_excel(visualization.activation_table_paths[0])
    second_table = pd.read_excel(visualization.activation_table_paths[1])
    assert list(first_table.columns) == ["Index", "X", "Y", "Z", "Activation"]
    assert first_table["Index"].tolist() == [1]
    assert first_table["Activation"].tolist() == [10]
    assert second_table["Activation"].tolist() == [-12]


def test_empty_activation_table_is_skipped_but_nifti_is_saved(tmp_path):
    result, coordinates = _visualizer_inputs()
    result.ActivationPatterns_BrainNetworks[:, 0] = 0

    visualization = BROADNESS_Visualizer(
        result,
        {
            "WhichPlots": [0, 0, 0, 0, 1],
            "MNI_coords": coordinates,
            "name_nii": tmp_path,
            "ncomps": [1],
        },
        show=False,
    )

    assert len(visualization.nifti_paths) == 1
    assert visualization.activation_table_paths == []
    nifti_values = np.asarray(nib.load(visualization.nifti_paths[0]).dataobj)
    assert np.count_nonzero(nifti_values) == 0


def test_nan_coordinate_rows_are_ignored_in_spatial_outputs():
    result, coordinates = _visualizer_inputs()
    coordinates[0] = np.nan

    visualization = BROADNESS_Visualizer(
        result,
        {
            "WhichPlots": [0, 0, 0, 1, 0],
            "MNI_coords": coordinates,
            "ncomps": [1],
        },
        show=False,
    )

    assert visualization.spatial_pattern_figure is not None


@pytest.mark.parametrize(
    ("options", "error", "message"),
    [
        ({"WhichPlots": [1, 1]}, ValueError, "five binary"),
        ({"WhichPlots": [0, 0, 2, 0, 0]}, ValueError, "zeros and ones"),
        (
            {"WhichPlots": [0, 0, 1, 0, 0], "ncomps": [0]},
            ValueError,
            "outside",
        ),
        (
            {"WhichPlots": [0, 0, 1, 0, 0], "ncomps": [1, 1]},
            ValueError,
            "duplicate",
        ),
        (
            {"WhichPlots": [0, 0, 1, 0, 0], "Labels": ["one"]},
            ValueError,
            "every condition",
        ),
        (
            {"WhichPlots": [0, 1, 0, 0, 0], "ncomps_var": 20},
            ValueError,
            "exceeds",
        ),
    ],
)
def test_invalid_options_are_rejected(options, error, message):
    result, _ = _visualizer_inputs()

    with pytest.raises(error, match=message):
        BROADNESS_Visualizer(result, options, show=False)


def test_spatial_outputs_require_matching_mni_coordinates():
    result, coordinates = _visualizer_inputs()

    with pytest.raises(ValueError, match="MNI_coords"):
        BROADNESS_Visualizer(
            result,
            {"WhichPlots": [0, 0, 0, 1, 0], "ncomps": [1]},
            show=False,
        )

    with pytest.raises(ValueError, match="sources x 3"):
        BROADNESS_Visualizer(
            result,
            {
                "WhichPlots": [0, 0, 0, 1, 0],
                "MNI_coords": coordinates[:-1],
                "ncomps": [1],
            },
            show=False,
        )


def test_nifti_output_requires_an_output_directory():
    result, coordinates = _visualizer_inputs()

    with pytest.raises(ValueError, match="name_nii"):
        BROADNESS_Visualizer(
            result,
            {
                "WhichPlots": [0, 0, 0, 0, 1],
                "MNI_coords": coordinates,
                "ncomps": [1],
            },
            show=False,
        )
