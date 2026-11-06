import numpy as np
import pytest
from scipy.stats import ttest_rel, f as f_distribution

from broadness import (
    BROADNESS_BehaviourCorrelation,
    BROADNESS_FDRCorrection,
    BROADNESS_RepeatedMeasuresANOVA,
    BROADNESS_SignificantSamplesToWindows,
)


def test_fdr_uses_one_family_and_retains_nan_positions():
    p = np.array([[0.001, np.nan, 0.02], [0.03, 0.2, np.nan]])
    significant, adjusted, critical = BROADNESS_FDRCorrection(p)
    np.testing.assert_array_equal(
        significant, [[True, False, True], [True, False, False]]
    )
    np.testing.assert_allclose(adjusted[~np.isnan(p)], [0.004, 0.04, 0.04, 0.2])
    assert np.isnan(adjusted[np.isnan(p)]).all()
    assert critical == 0.03


def test_fdr_empty_family_and_invalid_p_values():
    significant, adjusted, critical = BROADNESS_FDRCorrection([np.nan, np.nan])
    assert not significant.any()
    assert np.isnan(adjusted).all()
    assert np.isnan(critical)
    with pytest.raises(ValueError, match="between 0 and 1"):
        BROADNESS_FDRCorrection([0.1, 1.2])


def test_behaviour_correlation_uses_pairwise_valid_participants():
    behaviour = np.arange(1.0, 7.0)
    data = np.stack(
        [behaviour, behaviour[::-1], np.full(6, 3.0), [1, 2, np.nan, 4, 5, 6]]
    )
    r, p, adjusted, significant, critical = BROADNESS_BehaviourCorrelation(
        data, behaviour, correlationType="Spearman"
    )
    np.testing.assert_allclose(r[[0, 1, 3]], [1, -1, 1])
    assert np.isnan(r[2]) and np.isnan(p[2]) and np.isnan(adjusted[2])
    assert significant[[0, 1, 3]].all()
    assert not significant[2]
    assert critical <= 0.05


def test_two_condition_anova_matches_squared_paired_t():
    differences = np.array([1.0, 2.0, 1.5, 0.5, 2.5, 1.2])
    data = np.zeros((2, 1, 2, differences.size))
    data[0, 0, 0] = differences
    data[1, 0, 0] = np.array([1.0, -1.0, 0.5, -0.5, 0.2, -0.2])
    result = BROADNESS_RepeatedMeasuresANOVA(data)
    t_stat = ttest_rel(data[0, 0, 0], data[0, 0, 1]).statistic
    np.testing.assert_allclose(result.FStatistics[0, 0], t_stat**2)
    np.testing.assert_allclose(
        result.PValues[0, 0], f_distribution.sf(t_stat**2, 1, differences.size - 1)
    )
    assert result.GreenhouseGeisserEpsilon[0, 0] == 1
    assert result.Significant.shape == (2, 1)


def test_anova_drops_incomplete_participant_rows_and_applies_gg():
    generator = np.random.default_rng(11)
    data = generator.normal(size=(1, 1, 3, 8))
    data[0, 0, 2, 0] = np.nan
    result = BROADNESS_RepeatedMeasuresANOVA(data)
    assert 0.5 <= result.GreenhouseGeisserEpsilon[0, 0] <= 1
    np.testing.assert_allclose(
        result.ErrorDF[0, 0], result.GreenhouseGeisserEpsilon[0, 0] * 12
    )


def test_significant_samples_convert_to_contiguous_windows():
    time = np.linspace(-0.2, 0.6, 9)
    mask = [0, 1, 1, 0, 1, 0, 1, 1, 1]
    windows = BROADNESS_SignificantSamplesToWindows(mask, time, minimumSamples=2)
    assert windows == [(float(time[1]), float(time[2])), (float(time[6]), float(time[8]))]
    assert BROADNESS_SignificantSamplesToWindows(np.asarray(mask)[:, None], time[None, :], 2) == windows
