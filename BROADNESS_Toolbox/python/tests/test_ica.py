import numpy as np
import pytest

from broadness import BROADNESS_AlternativeNetworkEstimation_ICA, BROADNESSICAResult


def _mixed_sources(samples=800):
    generator = np.random.default_rng(12)
    sources = generator.laplace(size=(3, samples))
    mixing = np.array([[1.0, 0.3, -0.1], [0.2, 1.0, 0.4], [-0.3, 0.1, 1.0]])
    return mixing @ sources, sources, np.linspace(-0.2, 0.8, samples)


def test_jade_recovers_independent_sources_and_matlab_fields():
    data, sources, time = _mixed_sources()
    result = BROADNESS_AlternativeNetworkEstimation_ICA(data, time, icacomps=3)

    assert isinstance(result, BROADNESSICAResult)
    assert result.numICAcomps == 3
    assert result.TimeSeries_BrainNetworks.shape == (time.size, 3)
    assert result.ActivationPatterns_BrainNetworks.shape == (3, 3)
    np.testing.assert_array_equal(result.Time, time)
    np.testing.assert_array_equal(result.OriginalData, data)
    correlations = np.abs(np.corrcoef(result.TimeSeries_BrainNetworks.T, sources)[:3, 3:])
    assert np.all(np.max(correlations, axis=0) > 0.95)


def test_common_unmixing_is_applied_to_all_conditions_and_participants():
    data, _, time = _mixed_sources(samples=250)
    condition_offsets = np.array([0.0, 0.1])[None, None, :, None]
    participant_offsets = np.array([-0.05, 0.0, 0.08])[None, None, None, :]
    data = data[:, :, None, None] + condition_offsets + participant_offsets
    result = BROADNESS_AlternativeNetworkEstimation_ICA(data, time, icacomps=2)

    assert result.TimeSeries_BrainNetworks.shape == (250, 2, 2, 3)
    assert result.ActivationPatterns_BrainNetworks.shape == (3, 2)
    averaged = data.mean(axis=(2, 3))
    mean_scores = result.TimeSeries_BrainNetworks.mean(axis=(2, 3))
    separating = np.linalg.lstsq(averaged.T, mean_scores, rcond=None)[0].T
    for participant in range(3):
        for condition in range(2):
            np.testing.assert_allclose(
                result.TimeSeries_BrainNetworks[:, :, condition, participant],
                (separating @ data[:, :, condition, participant]).T,
                atol=1e-10,
            )
    np.testing.assert_allclose(
        result.ActivationPatterns_BrainNetworks,
        (separating @ np.cov(averaged)).T,
        atol=1e-10,
    )


def test_variance_threshold_selects_components_and_explicit_count_overrides_it():
    data, _, time = _mixed_sources()
    covariance_eigenvalues = np.linalg.eigvalsh(np.cov(data))[::-1]
    first_component_percent = 100 * covariance_eigenvalues[0] / covariance_eigenvalues.sum()
    result = BROADNESS_AlternativeNetworkEstimation_ICA(
        data, time, total_varexp=first_component_percent - 0.01
    )
    assert result.numICAcomps == 1

    explicit = BROADNESS_AlternativeNetworkEstimation_ICA(
        data, time, icacomps=2, total_varexp=1
    )
    assert explicit.numICAcomps == 2


@pytest.mark.parametrize(
    ("options", "error", "message"),
    [
        ({"icacomps": 0}, ValueError, "icacomps"),
        ({"icacomps": 4}, ValueError, "icacomps"),
        ({"icacomps": 1.5}, TypeError, "icacomps"),
        ({"total_varexp": 0}, ValueError, "total_varexp"),
        ({"total_varexp": 101}, ValueError, "total_varexp"),
        ({"total_varexp": "95"}, TypeError, "total_varexp"),
    ],
)
def test_invalid_options_are_rejected(options, error, message):
    data, _, time = _mixed_sources(samples=100)
    with pytest.raises(error, match=message):
        BROADNESS_AlternativeNetworkEstimation_ICA(data, time, **options)


def test_rank_deficiency_is_reported_before_jade_whitening():
    _, sources, time = _mixed_sources(samples=100)
    data = np.vstack([sources[0], sources[0] * 2, sources[1]])
    with pytest.raises(ValueError, match="rank"):
        BROADNESS_AlternativeNetworkEstimation_ICA(data, time, icacomps=3)
