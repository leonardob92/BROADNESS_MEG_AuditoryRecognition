import numpy as np
import pytest

from broadness import BROADNESSResult, BROADNESS_NetworkEstimation
from broadness.BROADNESS_NetworkEstimation import _randomize_data


def _synthetic_data():
    time = np.linspace(-0.2, 0.8, 31)
    source = np.arange(1, 7, dtype=float)[:, np.newaxis, np.newaxis, np.newaxis]
    sample = np.arange(time.size, dtype=float)[np.newaxis, :, np.newaxis, np.newaxis]
    condition = np.arange(1, 4, dtype=float)[np.newaxis, np.newaxis, :, np.newaxis]
    participant = np.arange(1, 5, dtype=float)[np.newaxis, np.newaxis, np.newaxis, :]

    data = (
        np.sin(0.17 * source * sample)
        + 0.2 * np.cos(0.11 * (source + condition) * sample)
        + 0.04 * source * condition
        + 0.03 * participant * np.sin(0.07 * sample)
        + 0.01 * participant * condition
    )
    return data, time


def test_two_dimensional_output_and_projection():
    data, time = _synthetic_data()
    data = data[:, :, 0, 0]

    result = BROADNESS_NetworkEstimation(data, time)

    assert isinstance(result, BROADNESSResult)
    assert result.ActivationPatterns_BrainNetworks.shape == (6, 6)
    assert result.Variance_BrainNetworks.shape == (6,)
    assert result.TimeSeries_BrainNetworks.shape == (31, 6)
    np.testing.assert_allclose(
        result.TimeSeries_BrainNetworks,
        data.T @ result.ActivationPatterns_BrainNetworks,
    )
    assert result.Significant_BrainNetworks == "MCS has not been run"
    assert result.VariancePermutations is None


def test_three_dimensional_conditions_are_preserved():
    data, time = _synthetic_data()
    data = data[:, :, :, 0]

    result = BROADNESS_NetworkEstimation(data, time)

    assert result.TimeSeries_BrainNetworks.shape == (31, 6, 3)
    for condition in range(3):
        np.testing.assert_allclose(
            result.TimeSeries_BrainNetworks[:, :, condition],
            data[:, :, condition].T @ result.ActivationPatterns_BrainNetworks,
        )


def test_four_dimensional_participants_use_common_patterns():
    data, time = _synthetic_data()

    result = BROADNESS_NetworkEstimation(data, time)
    direct_group_result = BROADNESS_NetworkEstimation(
        np.mean(data, axis=(2, 3)),
        time,
    )

    assert result.TimeSeries_BrainNetworks.shape == (31, 6, 3, 4)
    assert result.OriginalData.shape == data.shape
    np.testing.assert_allclose(
        result.ActivationPatterns_BrainNetworks,
        direct_group_result.ActivationPatterns_BrainNetworks,
    )
    for participant in range(4):
        for condition in range(3):
            np.testing.assert_allclose(
                result.TimeSeries_BrainNetworks[:, :, condition, participant],
                data[:, :, condition, participant].T
                @ result.ActivationPatterns_BrainNetworks,
            )


def test_time_window_uses_nearest_samples_and_includes_endpoints():
    data, time = _synthetic_data()

    result = BROADNESS_NetworkEstimation(
        data,
        time,
        time_window=[-0.09, 0.43],
    )

    start = int(np.argmin(np.abs(time - (-0.09))))
    end = int(np.argmin(np.abs(time - 0.43)))
    np.testing.assert_array_equal(result.Time, time[start : end + 1])
    np.testing.assert_array_equal(
        result.OriginalData,
        data[:, start : end + 1, :, :],
    )


@pytest.mark.parametrize("method", ["occurrences", "max_abs", "average"])
def test_sign_normalization_is_deterministic(method):
    data, time = _synthetic_data()

    first = BROADNESS_NetworkEstimation(data, time, sign_eigenvect=method)
    second = BROADNESS_NetworkEstimation(data, time, sign_eigenvect=method.upper())
    np.testing.assert_array_equal(
        first.ActivationPatterns_BrainNetworks,
        second.ActivationPatterns_BrainNetworks,
    )

    patterns = first.ActivationPatterns_BrainNetworks
    if method == "occurrences":
        assert np.all(np.mean(patterns > 0, axis=0) >= 0.5)
    elif method == "max_abs":
        indices = np.argmax(np.abs(patterns), axis=0)
        assert np.all(patterns[indices, np.arange(patterns.shape[1])] >= 0)
    else:
        assert np.all(np.mean(patterns, axis=0) >= -1e-15)


@pytest.mark.parametrize("randomization", [1, 2, 3])
def test_monte_carlo_is_reproducible_and_one_based(randomization):
    data, time = _synthetic_data()

    first = BROADNESS_NetworkEstimation(
        data,
        time,
        permutations_num=4,
        randomization=randomization,
        random_state=27,
    )
    second = BROADNESS_NetworkEstimation(
        data,
        time,
        permutations_num=4,
        randomization=randomization,
        random_state=27,
    )

    assert first.VariancePermutations is not None
    np.testing.assert_array_equal(
        first.VariancePermutations,
        second.VariancePermutations,
    )
    assert isinstance(first.Significant_BrainNetworks, np.ndarray)
    assert np.all(first.Significant_BrainNetworks >= 1)
    if first.Significant_BrainNetworks.size:
        np.testing.assert_array_equal(
            first.Significant_BrainNetworks,
            np.arange(1, first.Significant_BrainNetworks.size + 1),
        )
        assert first.TimeSeries_BrainNetworks.shape[1] == (
            first.Significant_BrainNetworks.size
        )
    else:
        assert first.TimeSeries_BrainNetworks.size == 0


@pytest.mark.parametrize("randomization", [1, 2, 3])
def test_randomization_preserves_the_expected_values(randomization):
    values = np.arange(35, dtype=float).reshape(5, 7)
    randomized = _randomize_data(
        values,
        randomization,
        np.random.default_rng(19),
    )

    assert randomized.shape == values.shape
    if randomization == 1:
        np.testing.assert_array_equal(
            np.sort(randomized, axis=1),
            np.sort(values, axis=1),
        )
    elif randomization == 2:
        np.testing.assert_array_equal(
            np.sort(randomized, axis=0),
            np.sort(values, axis=0),
        )
    else:
        np.testing.assert_array_equal(
            np.sort(randomized, axis=None),
            np.sort(values, axis=None),
        )


def test_negative_permutations_warn_and_skip_monte_carlo():
    data, time = _synthetic_data()

    with pytest.warns(UserWarning, match="negative"):
        result = BROADNESS_NetworkEstimation(data, time, permutations_num=-1)

    assert result.Significant_BrainNetworks == "MCS has not been run"
    assert result.VariancePermutations is None


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        ({"sign_eigenvect": "wrong"}, ValueError, "sign_eigenvect"),
        ({"sign_eigenvect": 1}, TypeError, "sign_eigenvect"),
        ({"randomization": 0}, ValueError, "randomization"),
        ({"randomization": 4}, ValueError, "randomization"),
        ({"randomization": 1.5}, TypeError, "randomization"),
        ({"permutations_num": 1.5}, TypeError, "permutations_num"),
        ({"random_state": 2.5}, TypeError, "random_state"),
        ({"random_state": -1}, ValueError, "random_state"),
        ({"time_window": [0.5, 0.1]}, ValueError, "time_window"),
        ({"time_window": [0.1]}, ValueError, "time_window"),
    ],
)
def test_invalid_options_are_rejected(kwargs, error, message):
    data, time = _synthetic_data()

    with pytest.raises(error, match=message):
        BROADNESS_NetworkEstimation(data, time, **kwargs)


def test_invalid_data_and_time_are_rejected():
    data, time = _synthetic_data()

    with pytest.raises(ValueError, match="shape"):
        BROADNESS_NetworkEstimation(data[np.newaxis, ...], time)
    with pytest.raises(ValueError, match="equal the length"):
        BROADNESS_NetworkEstimation(data, time[:-1])

    invalid = data.copy()
    invalid[0, 0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        BROADNESS_NetworkEstimation(invalid, time)


def test_zero_variance_data_are_rejected():
    time = np.linspace(0, 1, 10)
    data = np.ones((4, 10))

    with pytest.raises(ValueError, match="zero variance"):
        BROADNESS_NetworkEstimation(data, time)
