"""Small headless checks of figure modes and exported files."""

from types import SimpleNamespace
import pickle

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
import numpy as np
import pytest

from broadness import (
    BROADNESS_FigureSettings,
    BROADNESS_FinalizeFigure,
    BROADNESS_PhaseSpace_RQA,
    BROADNESS_Plot_ActivationTimeseries,
    BROADNESS_SpatialActivationClustering,
    BROADNESS_SpatialPatternGroupStatistics,
    BROADNESS_SpatialPatternStatistics,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _data():
    generator = np.random.default_rng(8)
    coords = np.array([(x, y, z) for x in (0, 8) for y in (0, 8) for z in (0, 8)], dtype=float)
    maps = generator.normal(size=(8, 2, 2, 6))
    scores = generator.normal(size=(7, 2, 2, 6))
    return SimpleNamespace(
        OriginalData=generator.normal(size=(8, 7, 2, 6)),
        ActivationPatterns_BrainNetworks=maps.mean(axis=(2, 3)),
        ActivationPatterns_BrainNetworks_Individual=maps,
        TimeSeries_BrainNetworks=scores, Time=np.linspace(-0.2, 0.4, 7),
    ), coords


def test_settings_and_finalize_export(tmp_path):
    settings = BROADNESS_FigureSettings("save", "both", tmp_path, ("png", "pdf", "fig"), "Study", "Test")
    assert settings.MakeIndividual and settings.MakeSummary
    fig, ax = plt.subplots()
    ax.plot([0, 1], [1, 0])
    files = BROADNESS_FinalizeFigure(fig, settings, "Example")
    assert {path.suffix for path in files} == {".png", ".pdf", ".fig"}
    assert all(path.exists() and path.stat().st_size > 100 for path in files)
    with next(path for path in files if path.suffix == ".fig").open("rb") as stream:
        loaded = pickle.load(stream)
    assert loaded.axes[0].lines[0].get_ydata().tolist() == [1, 0]
    plt.close(loaded)
    assert not plt.fignum_exists(fig.number)
    with pytest.raises(ValueError, match="outputPath"):
        BROADNESS_FigureSettings("save")


def test_activation_timeseries_individual_and_summary(tmp_path):
    broadness, _ = _data()
    figures = BROADNESS_Plot_ActivationTimeseries(
        broadness.TimeSeries_BrainNetworks, broadness.Time,
        Groups=([1, 2, 3], [4, 5, 6]),
        ConditionLabels=["Old", "New"], GroupLabels=["A", "B"],
        SignificantWindows=[[(0, 0.1)], [(0.2, 0.3)]],
        SignificanceLabels=["contrast", "interaction"],
        SignificanceStyle="line", STEStyle=2,
        FigureMode="save", FigureLayout="both", OutputPath=tmp_path,
    )
    assert figures.Handles == []
    assert len(figures.Files) == 3
    assert all(path.exists() for path in figures.Files)


def test_analysis_figure_modes_save_small_outputs(tmp_path):
    broadness, coords = _data()
    rqa, rqa_figures = BROADNESS_PhaseSpace_RQA(
        broadness, principalcomps=(1, 2), threshold=0.2,
        figuremode="save", figurelayout="summary", OutputPath=tmp_path,
    )
    assert len(rqa_figures.Files) == 2
    assert rqa.Figures is rqa_figures
    paired, paired_figures = BROADNESS_SpatialPatternStatistics(
        broadness, principalcomps=(1,), mni_coords=coords, permutations=3,
        figuremode="save", figurelayout="summary", OutputPath=tmp_path,
    )
    assert len(paired_figures.Files) == 1
    assert paired.Figures is paired_figures
    groups, group_figures = BROADNESS_SpatialPatternGroupStatistics(
        broadness, [1, 1, 1, 2, 2, 2], principalcomps=(1,),
        mni_coords=coords, permutations=3,
        figuremode="save", figurelayout="summary", OutputPath=tmp_path,
    )
    assert len(group_figures.Files) == 2
    assert groups.Figures is group_figures
    cluster, clustering_figures = BROADNESS_SpatialActivationClustering(
        broadness, principalcomps=(1, 2), nclusters=(2,),
        mni_coords=coords, kmeans_restarts=2,
        figuremode="save", figurelayout="summary", figureoutpath=tmp_path,
    )
    assert len(clustering_figures.Files) == 1
    assert cluster.Figures is clustering_figures
    for files in (rqa_figures.Files, paired_figures.Files, group_figures.Files, clustering_figures.Files):
        assert all(path.exists() for path in files)


def test_analysis_individual_figure_branches(tmp_path):
    broadness, coords = _data()
    broadness.TimeSeries_BrainNetworks = broadness.TimeSeries_BrainNetworks[:, :, :1]
    _, rqa_figures = BROADNESS_PhaseSpace_RQA(
        broadness, figuremode="save", figurelayout="individual", OutputPath=tmp_path,
    )
    assert len(rqa_figures.Files) == 3
    broadness.ActivationPatterns_BrainNetworks_Individual = broadness.ActivationPatterns_BrainNetworks_Individual[:, :, :1]
    _, group_figures = BROADNESS_SpatialPatternGroupStatistics(
        broadness, [1, 1, 1, 2, 2, 2], principalcomps=(1,),
        mni_coords=coords, permutations=3,
        figuremode="save", figurelayout="individual", OutputPath=tmp_path,
    )
    assert len(group_figures.Files) == 1
    _, clustering_figures = BROADNESS_SpatialActivationClustering(
        broadness, nclusters=(2,), mni_coords=coords, kmeans_restarts=2,
        figuremode="save", figurelayout="individual", figureoutpath=tmp_path,
    )
    assert len(clustering_figures.Files) == 4


def test_phase_video_returns_live_animation():
    broadness, _ = _data()
    broadness.TimeSeries_BrainNetworks = broadness.TimeSeries_BrainNetworks[:, :, :1]
    with pytest.warns(UserWarning, match="non-interactive"):
        result, figures = BROADNESS_PhaseSpace_RQA(broadness, video="on")
    assert result.Figures is figures
    assert len(figures.Animations) == 1
    assert len(figures.Handles) == 1
    figures.Animations[0]._draw_next_frame(0, blit=False)
