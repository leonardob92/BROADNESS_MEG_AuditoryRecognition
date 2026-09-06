import numpy as np
import pytest

from broadness import BROADNESS_EffectiveDimensionality


def test_equal_components_return_their_number():
    assert BROADNESS_EffectiveDimensionality([25, 25, 25, 25]) == 4


def test_single_dominant_component_returns_one():
    assert BROADNESS_EffectiveDimensionality([100, 0, 0, 0]) == 1


def test_typical_eigenspectrum_returns_participation_ratio():
    assert BROADNESS_EffectiveDimensionality([50, 30, 20]) == 3


def test_multiple_eigenspectra_are_computed_along_first_dimension():
    eigenspectra = np.array(
        [
            [40, 100, 30],
            [30, 0, 30],
            [20, 0, 30],
            [10, 0, 30],
        ],
        dtype=float,
    )

    result = BROADNESS_EffectiveDimensionality(eigenspectra)

    np.testing.assert_array_equal(result, np.array([3, 1, 4]))


def test_singleton_dimensions_are_removed_like_matlab_squeeze():
    eigenspectra = np.full((4, 1, 2), 25.0)

    result = BROADNESS_EffectiveDimensionality(eigenspectra)

    np.testing.assert_array_equal(result, np.array([4, 4]))


def test_result_is_invariant_to_positive_scaling():
    eigenspectrum = np.array([60, 25, 10, 5], dtype=float)

    original = BROADNESS_EffectiveDimensionality(eigenspectrum)
    scaled = BROADNESS_EffectiveDimensionality(eigenspectrum * 1000)

    assert original == scaled


def test_zero_eigenspectrum_uses_matlab_realmin_guard():
    assert BROADNESS_EffectiveDimensionality([0, 0, 0]) == 0


def test_positive_half_integer_uses_matlab_rounding():
    eigenspectrum = np.array([1.0, 2.0 - np.sqrt(3.0)])

    assert BROADNESS_EffectiveDimensionality(eigenspectrum) == 2


@pytest.mark.parametrize(
    ("eigenspectrum", "error", "message"),
    [
        (1.0, ValueError, "at least one component"),
        ([], ValueError, "empty"),
        ([1, -1], ValueError, "negative"),
        ([1, np.nan], ValueError, "NaN"),
        ([1, np.inf], ValueError, "infinite"),
        ([1, 2j], TypeError, "real numeric"),
        (["one", "two"], TypeError, "real numeric"),
    ],
)
def test_invalid_eigenspectra_are_rejected(eigenspectrum, error, message):
    with pytest.raises(error, match=message):
        BROADNESS_EffectiveDimensionality(eigenspectrum)
