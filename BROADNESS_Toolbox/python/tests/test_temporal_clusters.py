import numpy as np
import pandas as pd
import pytest
from scipy.stats import ttest_ind

from broadness import (
    BROADNESS_TemporalClusterCorrection,
    BROADNESS_TemporalClusterCorrection_ANOVA,
    BROADNESS_TemporalClusterCorrection_GroupANOVA,
    BROADNESS_TemporalClusterCorrection_Groups,
)
from broadness._temporal_clusters import _clusters, _welch_t


def _effect_data():
    generator = np.random.default_rng(14)
    values = generator.normal(0, 0.4, size=(20, 1, 2, 24))
    values[6:12, 0, 0, :8] += 3.0
    return values, np.repeat([1, 2, 3], 8), np.arange(20) * 0.01


def test_cluster_labels_keep_signs_separate_and_tables_are_chronological():
    labels, masses = _clusters(np.array([-5.0, -6.0, 0.0, 7.0, 8.0]), 4.0, two_sided=True)
    np.testing.assert_array_equal(labels, [2, 2, 0, 1, 1])
    np.testing.assert_array_equal(masses, [15.0, 11.0])


def test_paired_sign_flips_find_effect_and_report_cluster_level_p():
    data, _, times = _effect_data()
    first = BROADNESS_TemporalClusterCorrection(
        data, [[1, 2]], permutations=39, randomseed=2, time=times
    )
    second = BROADNESS_TemporalClusterCorrection(
        data, [[1, 2]], permutations=39, randomseed=2, time=times
    )
    significant, corrected, critical, details = first
    assert significant.shape == (20, 1, 1)
    assert corrected.shape == significant.shape
    assert details.NullMaximumClusterMass.shape == (39, 1, 1)
    np.testing.assert_array_equal(
        details.NullMaximumClusterMass, second[3].NullMaximumClusterMass
    )
    assert significant[6:12, 0, 0].all()
    assert critical == 0.05
    assert isinstance(details.Clusters[0, 0], pd.DataFrame)
    assert list(details.Clusters[0, 0].columns[:3]) == ["Cluster", "StartTime", "EndTime"]
    assert details.Clusters[0, 0].iloc[0].StartTime == times[6]
    assert np.unique(corrected[6:12, 0, 0]).size == 1
    assert np.isnan(corrected[0, 0, 0])
    assert "separately per network and condition pair" in details.Settings.CorrectionFamily


def test_anova_permutation_uses_one_null_per_network():
    data, _, _ = _effect_data()
    significant, corrected, critical, details = BROADNESS_TemporalClusterCorrection_ANOVA(
        data, permutations=39, randomseed=2
    )
    assert significant.shape == (20, 1)
    assert details.NullMaximumClusterMass.shape == (39, 1)
    assert significant[6:12, 0].all()
    assert critical == 0.05
    np.testing.assert_array_equal(details.DegreesFreedom, [1, 23])
    assert list(details.Clusters[0].columns[:3]) == ["Cluster", "StartSample", "EndSample"]
    assert details.Clusters[0].iloc[0].StartSample == 7
    assert details.FStatistics[8, 0] > details.ClusterFormingThreshold
    assert np.isnan(corrected[0, 0])


def test_group_wild_bootstrap_shapes_and_detects_planted_effect():
    data, labels, times = _effect_data()
    significant, corrected, _, details = BROADNESS_TemporalClusterCorrection_Groups(
        data, labels, [[1, 2]], permutations=39, randomseed=2, time=times
    )
    assert significant.shape == (20, 1, 2, 1)
    assert details.ClusterFormingThreshold.shape == significant.shape
    assert details.NullMaximumClusterMass.shape == (39, 1, 2, 1)
    assert significant[6:12, 0, 0, 0].all()
    assert details.Clusters[0, 0, 0].iloc[0].StartTime == times[6]
    assert corrected[7, 0, 0, 0] <= 0.05
    assert "wild bootstrap" in details.Settings.Resampling


def test_group_anova_uses_welch_f_and_separate_condition_families():
    data, labels, _ = _effect_data()
    significant, corrected, _, details = BROADNESS_TemporalClusterCorrection_GroupANOVA(
        data, labels, permutations=39, randomseed=2
    )
    assert significant.shape == (20, 1, 2)
    assert details.NullMaximumClusterMass.shape == (39, 1, 2)
    np.testing.assert_array_equal(details.Settings.GroupIDs, [1, 2, 3])
    assert significant[6:12, 0, 0].all()
    assert details.Clusters[0, 0].iloc[0].StartSample == 7
    assert corrected[7, 0, 0] <= 0.05


def test_welch_t_matches_scipy_for_each_time_sample():
    rng = np.random.default_rng(9)
    x = rng.normal(size=(6, 4))
    y = rng.normal(loc=0.4, scale=1.5, size=(6, 7))
    t_values, df = _welch_t(x, y)
    reference = ttest_ind(x, y, axis=1, equal_var=False)
    np.testing.assert_allclose(t_values, reference.statistic)
    np.testing.assert_allclose(df, reference.df)


def test_invalid_pairs_and_group_sizes_are_rejected():
    data, labels, _ = _effect_data()
    with pytest.raises(ValueError, match="condition_pairs"):
        BROADNESS_TemporalClusterCorrection(data, [[1, 2], [2, 1]], permutations=2)
    with pytest.raises(ValueError, match="group_pairs"):
        BROADNESS_TemporalClusterCorrection_Groups(data, labels, [[1, 4]], permutations=2)
    with pytest.raises(ValueError, match="three groups"):
        BROADNESS_TemporalClusterCorrection_GroupANOVA(
            data, np.repeat([1, 2], 12), permutations=2
        )
    with pytest.raises(ValueError, match="permutations"):
        BROADNESS_TemporalClusterCorrection_ANOVA(data, permutations=0)
