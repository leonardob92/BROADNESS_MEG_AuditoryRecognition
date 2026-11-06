"""Small numerical fixtures for the phase-space and spatial analysis ports."""

from types import SimpleNamespace

import numpy as np
import pytest

from broadness import (
    BROADNESS_PhaseSpace_RQA,
    BROADNESS_PhaseSpaceStatistics,
    BROADNESS_SpatialActivationClustering,
    BROADNESS_SpatialPatternGroupStatistics,
    BROADNESS_SpatialPatternStatistics,
)


def _fixture():
    generator = np.random.default_rng(14)
    time = np.linspace(-0.2, 0.5, 12)
    scores = generator.normal(size=(12, 3, 3, 9))
    scores[:, 0, 1] += np.linspace(0, 1, 12)[:, None]
    scores[:, 1, 2] -= 0.5
    coords = np.array([(x, y, z) for x in (0, 8) for y in (0, 8) for z in (0, 8)], dtype=float)
    maps = generator.normal(size=(8, 3, 3, 9))
    maps[:4, 0, 1] += 1.5
    maps[4:, 1, 2, 6:] -= 1
    mean_maps = maps.mean(axis=(2, 3))
    return SimpleNamespace(
        Time=time, TimeSeries_BrainNetworks=scores,
        ActivationPatterns_BrainNetworks=mean_maps,
        ActivationPatterns_BrainNetworks_Individual=maps,
        OriginalData=generator.normal(size=(8, 12, 3, 9)),
    ), coords


def test_phase_rqa_and_statistics():
    broadness, _ = _fixture()
    rqa, figures = BROADNESS_PhaseSpace_RQA(
        broadness, principalcomps=(1, 2), timeinterval=(-0.1, 0.4),
        threshold=0.3, normalization="pooled_zscore", theiler_window=1,
    )
    assert figures.Files == []
    assert rqa.PhaseSpace.ParticipantCoordinates.shape == (3, 9)
    assert len(rqa.RQA_metrics) == 9
    assert list(rqa.RQA_metrics[0].columns) == ["RR", "L", "DET", "ENTR", "TT", "LAM", "V_max", "DIV"]
    distance = rqa.RecurrencePlots.DistMat[0, 0]
    recurrence = rqa.RecurrencePlots.RecurPlot[0, 0]
    np.testing.assert_allclose(distance, distance.T)
    assert not np.any(np.diag(recurrence))
    assert not np.any(np.diag(recurrence, 1))
    assert np.isclose(rqa.RQA_metrics[0].RR.iloc[0], recurrence.sum() / (recurrence.size - 3 * len(recurrence) + 2))

    stats = BROADNESS_PhaseSpaceStatistics(
        rqa, groupindices=[list(range(1, 5)), list(range(5, 10))],
        testtype="both", testmetrics=("Speed", "TotalPathLength"),
    )
    assert stats.Trajectories.shape == (len(rqa.PhaseSpace.Time), 2, 3, 9)
    assert stats.Metrics.Speed.shape == (len(rqa.PhaseSpace.Time), 3, 9)
    assert stats.Tests.Enabled
    assert len(stats.Tests.PhaseSpace.ConditionComparisons) == 9
    assert len(stats.Tests.PhaseSpace.GroupComparisons) == 3
    assert stats.Tests.Metrics.Speed.Omnibus.Parametric.Condition.F.shape == (len(rqa.PhaseSpace.Time),)
    assert stats.Tests.Metrics.TotalPathLength.ConditionComparisons[0].Parametric.P.shape == (1,)


def test_rqa_full_recurrence_has_expected_line_statistics():
    scores = np.arange(4, dtype=float)[:, None]
    rqa, _ = BROADNESS_PhaseSpace_RQA(
        SimpleNamespace(Time=np.arange(4), TimeSeries_BrainNetworks=scores),
        principalcomps=(1,), threshold=1.1,
    )
    metrics = rqa.RQA_metrics[0].iloc[0]
    assert metrics.RR == 1
    assert np.isclose(metrics.L, (2 + 2 + 3 + 3 + 4) / 5)
    assert np.isclose(metrics.DET, 14 / 16)
    assert metrics.TT == 4
    assert metrics.LAM == 1
    assert metrics.V_max == 4
    assert np.isclose(metrics.DIV, 1 / 3)


def test_spatial_paired_and_group():
    broadness, coords = _fixture()
    paired, figures = BROADNESS_SpatialPatternStatistics(
        broadness, principalcomps=(1, 2), mni_coords=coords,
        omnibus="on", permutations=11, clusteralpha=0.2, randomseed=3,
    )
    assert figures.Files == []
    assert paired.Voxelwise.Pairwise.TStatistics.shape == (8, 2, 3)
    assert paired.Voxelwise.Pairwise.NullMaximumClusterMass.shape == (11, 3)
    assert paired.Voxelwise.Omnibus.FStatistics.shape == (8, 2)
    assert paired.Summary.Values.shape == (6, 2, 3, 9)
    assert paired.SummaryTests.Pairwise.AdjustedPValues.shape == (6, 2, 3)
    np.testing.assert_allclose(paired.Summary.MapStrength, np.sqrt(np.mean(broadness.ActivationPatterns_BrainNetworks_Individual[:, :2] ** 2, axis=0)))

    groups, _ = BROADNESS_SpatialPatternGroupStatistics(
        broadness, [1] * 3 + [2] * 3 + [3] * 3,
        principalcomps=(1, 2), mni_coords=coords,
        permutations=11, clusteralpha=0.2, randomseed=3,
    )
    assert groups.Voxelwise.Pairwise.TStatistics.shape == (8, 2, 3, 3)
    assert groups.Voxelwise.Omnibus.FStatistics.shape == (8, 2, 3)
    assert groups.Voxelwise.Pairwise.NullMaximumClusterMass.shape == (11, 3, 3)
    assert groups.Summary.Values.shape == (6, 2, 3, 9)

    two_groups, _ = BROADNESS_SpatialPatternGroupStatistics(
        broadness, [1] * 4 + [2] * 5,
        principalcomps=(1, 2), mni_coords=coords, permutations=3,
    )
    assert not two_groups.Voxelwise.Omnibus.Enabled
    assert two_groups.Voxelwise.Pairwise.GroupPairs.shape == (1, 2)


def test_spatial_clustering_result_structure():
    broadness, coords = _fixture()
    result, figures = BROADNESS_SpatialActivationClustering(
        broadness, principalcomps=(1, 2), nclusters=(2, 3),
        evalclusters=1, mni_coords=coords, kmeans_restarts=2,
    )
    assert figures.Files == []
    assert result.optimalK in (2, 3)
    assert result.idx.shape == (8, 2)
    assert len(result.Clusters_info) == result.optimalK
    assert result.ClusterSummary.VoxelCount.sum() == 8
    assert sorted(result.idx[f"Nclusters_{result.optimalK}"].unique()) == list(range(1, result.optimalK + 1))
    with pytest.raises(ValueError, match="sources in the bundled MNI voxel order"):
        BROADNESS_SpatialActivationClustering(
            broadness, principalcomps=(1, 2), nclusters=(2,),
            mni_coords=coords, kmeans_restarts=2, OutputPath="unused",
        )
