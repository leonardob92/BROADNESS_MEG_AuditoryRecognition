"""PCA-based broadband brain network estimation for BROADNESS.

=========================================================================
 USER DOCUMENTATION
=========================================================================

This module provides the public function ``BROADNESS_NetworkEstimation``.
It takes a multivariate M/EEG dataset and estimates broadband brain networks
using principal component analysis (PCA).

INPUT DATA
----------
``data`` can contain:

- ``sources x time``
- ``sources x time x conditions``
- ``sources x time x conditions x participants``

For participant-level data, BROADNESS first averages across participants and
conditions to estimate a common set of group-level spatial activation
patterns. Every participant's original data are then projected onto those
common patterns to obtain participant-specific network time series.

MAIN OPTIONS
------------
``time_window`` selects an optional interval in seconds.
``permutations_num`` requests Monte Carlo simulations.
``randomization`` selects time (1), space (2), or combined (3) randomization.
``sign_eigenvect`` selects ``"occurrences"``, ``"max_abs"``, or ``"average"``.
``random_state`` optionally makes the Python randomization reproducible.

EXAMPLE
-------
Run the default analysis with::

    from broadness import BROADNESS_NetworkEstimation

    BROADNESS = BROADNESS_NetworkEstimation(data, time)

Run it with optional settings using::

    BROADNESS = BROADNESS_NetworkEstimation(
        data,
        time,
        time_window=[0.350, 1.750],
        permutations_num=100,
        randomization=1,
        sign_eigenvect="max_abs",
        random_state=1,
    )

MAIN OUTPUTS
------------
``BROADNESS.ActivationPatterns_BrainNetworks`` contains spatial patterns.
``BROADNESS.Variance_BrainNetworks`` contains explained variance percentages.
``BROADNESS.TimeSeries_BrainNetworks`` contains the network time series.
``BROADNESS.Significant_BrainNetworks`` contains one-based network numbers
when Monte Carlo simulations are requested.
``BROADNESS.VariancePermutations`` contains the mean randomized eigenspectrum,
or ``None`` when Monte Carlo simulations are not requested.
``BROADNESS.Time`` and ``BROADNESS.OriginalData`` contain the analysed time
vector and input data within the selected time window.

INDEXING
--------
Public component and participant selections use one-based numbers for
consistency with MATLAB. Direct NumPy array indexing remains zero-based.

``BROADNESS_NetworkEstimation`` is the analysis function intended for users.
It returns a ``BROADNESSResult`` object containing the analysis outputs.
Functions whose names begin with an underscore are internal computational
helpers and do not need to be called by users.

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

from dataclasses import dataclass
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import svd


# =========================================================================
#  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS)
#  NETWORK ESTIMATION
# =========================================================================

FloatArray = NDArray[np.float64]
IntegerArray = NDArray[np.int64]


@dataclass(slots=True)
class BROADNESSResult:
    """Outputs of :func:`BROADNESS_NetworkEstimation`.

    Array dimensions intentionally mirror the MATLAB toolbox. Original data
    retain ``sources x time x [conditions] x [participants]`` dimensions,
    while network time series use
    ``time x components x [conditions] x [participants]``.
    """

    Variance_BrainNetworks: FloatArray
    Significant_BrainNetworks: IntegerArray | str
    ActivationPatterns_BrainNetworks: FloatArray
    TimeSeries_BrainNetworks: FloatArray
    Time: FloatArray
    OriginalData: NDArray[np.number]
    VariancePermutations: FloatArray | None = None


# =========================================================================
#  INTERNAL HELPER FUNCTIONS
#  Users do not need to call the underscore-prefixed functions below.
# =========================================================================


def _numeric_data(data: ArrayLike) -> NDArray[np.number]:
    """Validate the documented two-, three-, or four-dimensional input."""
    values = np.asarray(data)
    if values.ndim not in (2, 3, 4):
        raise ValueError(
            "data must have shape sources x time, sources x time x "
            "conditions, or sources x time x conditions x participants"
        )
    if values.shape[0] < 2 or values.shape[1] < 2:
        raise ValueError("data must contain at least two sources and two time points")
    if any(size == 0 for size in values.shape):
        raise ValueError("data dimensions cannot be empty")
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("data must contain real numeric values")
    if not np.all(np.isfinite(values)):
        raise ValueError("data contains NaN or infinite values")
    return values


def _time_vector(time: ArrayLike, n_time: int) -> FloatArray:
    """Validate and return time as a one-dimensional floating-point vector."""
    values = np.asarray(time)
    if values.ndim == 0 or values.ndim > 2:
        raise ValueError("time must be a vector")
    if values.ndim == 2 and 1 not in values.shape:
        raise ValueError("time must be a vector")
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("time must contain real numeric values")

    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size != n_time:
        raise ValueError(
            "the number of time points in data must equal the length of time"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("time contains NaN or infinite values")
    return values


def _time_indices(time: FloatArray, time_window: ArrayLike | None) -> tuple[int, int]:
    """Return inclusive indices for the requested nearest-sample time window."""
    if time_window is None:
        print("Duration not provided.. using full data duration")
        return 0, time.size - 1

    window = np.asarray(time_window)
    if (
        window.ndim != 1
        or window.size != 2
        or not np.issubdtype(window.dtype, np.number)
        or np.iscomplexobj(window)
    ):
        raise ValueError("time_window must contain two real numeric values")
    window = np.asarray(window, dtype=float)
    if not np.all(np.isfinite(window)):
        raise ValueError("time_window contains NaN or infinite values")

    index_start = int(np.argmin(np.abs(time - window[0])))
    index_end = int(np.argmin(np.abs(time - window[1])))
    if index_end <= index_start:
        raise ValueError("the end of time_window must be after its beginning")
    return index_start, index_end


def _positive_integer(value: int, name: str, *, allow_zero: bool) -> int:
    """Validate an integer option without silently truncating it."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    integer = int(value)
    minimum = 0 if allow_zero else 1
    if integer < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return integer


def _compute_pca(data: NDArray[np.number]) -> tuple[FloatArray, FloatArray]:
    """Compute MATLAB-compatible PCA coefficients and explained variance."""
    observations = np.asarray(data, dtype=float).T
    observations = observations - np.mean(observations, axis=0, keepdims=True)

    number_observations, number_variables = observations.shape
    degrees_freedom = number_observations - 1
    number_components = min(number_variables, degrees_freedom)
    if number_components < 1:
        raise ValueError("PCA requires at least two analysed time points")

    _, singular_values, right_vectors = svd(
        observations,
        full_matrices=False,
        check_finite=False,
        lapack_driver="gesdd",
    )
    singular_values = singular_values[:number_components]
    activation_patterns = right_vectors[:number_components, :].T
    eigenvalues = singular_values**2 / degrees_freedom

    total_variance = np.sum(eigenvalues)
    if not np.isfinite(total_variance) or total_variance <= 0:
        raise ValueError("PCA cannot be computed because the data have zero variance")
    variance = eigenvalues * 100.0 / total_variance
    return activation_patterns, variance


def _normalize_eigenvector_signs(
    activation_patterns: FloatArray,
    method: str,
) -> FloatArray:
    """Apply one of the three sign conventions used by MATLAB BROADNESS."""
    normalized = activation_patterns.copy()

    if method == "occurrences":
        signs = np.where(np.mean(normalized > 0, axis=0) < 0.5, -1.0, 1.0)
    elif method == "max_abs":
        maximum_indices = np.argmax(np.abs(normalized), axis=0)
        maximum_values = normalized[maximum_indices, np.arange(normalized.shape[1])]
        signs = np.where(maximum_values < 0, -1.0, 1.0)
    else:
        signs = np.where(np.mean(normalized, axis=0) < 0, -1.0, 1.0)

    normalized *= signs
    return normalized


def _randomize_data(
    averaged_data: NDArray[np.number],
    randomization: int,
    generator: np.random.Generator,
) -> FloatArray:
    """Randomize time, space, or both dimensions following MATLAB BROADNESS."""
    values = np.asarray(averaged_data, dtype=float)
    randomized = np.empty_like(values)

    if randomization == 1:
        for source in range(values.shape[0]):
            randomized[source, :] = values[source, generator.permutation(values.shape[1])]
    elif randomization == 2:
        for sample in range(values.shape[1]):
            randomized[:, sample] = values[generator.permutation(values.shape[0]), sample]
    else:
        flattened = values.ravel(order="F")
        randomized_flattened = flattened[generator.permutation(flattened.size)]
        randomized = randomized_flattened.reshape(values.shape, order="F")

    return randomized


def _network_time_series(
    original_data: NDArray[np.number],
    activation_patterns: FloatArray,
    number_components: int,
) -> FloatArray:
    """Project every condition and participant onto the common patterns."""
    selected_patterns = activation_patterns[:, :number_components]

    if original_data.ndim == 2:
        return np.asarray(original_data.T @ selected_patterns, dtype=float)

    if original_data.ndim == 3:
        time_series = np.empty(
            (original_data.shape[1], number_components, original_data.shape[2]),
            dtype=float,
        )
        for condition in range(original_data.shape[2]):
            time_series[:, :, condition] = (
                original_data[:, :, condition].T @ selected_patterns
            )
        return time_series

    time_series = np.empty(
        (
            original_data.shape[1],
            number_components,
            original_data.shape[2],
            original_data.shape[3],
        ),
        dtype=float,
    )
    for participant in range(original_data.shape[3]):
        print(
            "Computing time series for participant "
            f"{participant + 1} / {original_data.shape[3]}"
        )
        for condition in range(original_data.shape[2]):
            time_series[:, :, condition, participant] = (
                original_data[:, :, condition, participant].T @ selected_patterns
            )
    return time_series


# =========================================================================
#  MAIN PUBLIC FUNCTION
# =========================================================================


def BROADNESS_NetworkEstimation(
    data: ArrayLike,
    time: ArrayLike,
    *,
    time_window: ArrayLike | None = None,
    permutations_num: int = 0,
    randomization: int = 1,
    sign_eigenvect: str = "max_abs",
    random_state: int | None = None,
) -> BROADNESSResult:
    """Estimate broadband brain networks using PCA.

    Parameters
    ----------
    data
        Data with dimensions ``sources x time x [conditions] x
        [participants]``. For participant-level data, a common decomposition
        is estimated after averaging across participants and conditions. Each
        participant and condition is then projected onto those common spatial
        activation patterns.
    time
        Time vector in seconds, with one value per sample in ``data``.
    time_window
        Optional ``[start, end]`` interval in seconds. The closest available
        samples are selected and both endpoints are included.
    permutations_num
        Number of Monte Carlo randomizations. The default is zero.
    randomization
        Randomization strategy: 1 for time, 2 for space, or 3 for both.
    sign_eigenvect
        Eigenvector sign convention: ``"occurrences"``, ``"max_abs"``, or
        ``"average"``. The default is ``"max_abs"``.
    random_state
        Optional seed for reproducible Python Monte Carlo randomization.

    Returns
    -------
    BROADNESSResult
        Spatial activation patterns, explained variance, network time series,
        time, analysed original data, and optional Monte Carlo results.

    Notes
    -----
    Significant brain-network numbers are one-based for consistency with the
    MATLAB toolbox. Direct NumPy array indexing remains zero-based.
    """
    print("Checking inputs")

    values = _numeric_data(data)
    time_values = _time_vector(time, values.shape[1])

    if not isinstance(sign_eigenvect, str):
        raise TypeError(
            'sign_eigenvect must be "occurrences", "max_abs", or "average"'
        )
    sign_method = sign_eigenvect.lower()
    valid_sign_options = {"occurrences", "max_abs", "average"}
    if sign_method not in valid_sign_options:
        raise ValueError(
            'sign_eigenvect must be "occurrences", "max_abs", or "average"'
        )

    randomization = _positive_integer(randomization, "randomization", allow_zero=False)
    if randomization not in (1, 2, 3):
        raise ValueError("randomization must be 1, 2, or 3")

    if isinstance(permutations_num, (int, np.integer)) and int(permutations_num) < 0:
        warnings.warn(
            "permutations_num is negative; Monte Carlo simulations will not be run",
            UserWarning,
            stacklevel=2,
        )
        permutations = 0
    else:
        permutations = _positive_integer(
            permutations_num,
            "permutations_num",
            allow_zero=True,
        )

    if random_state is not None and (
        isinstance(random_state, (bool, np.bool_))
        or not isinstance(random_state, (int, np.integer))
    ):
        raise TypeError("random_state must be an integer or None")
    if random_state is not None and int(random_state) < 0:
        raise ValueError("random_state must be non-negative")

    index_start, index_end = _time_indices(time_values, time_window)
    analysed_data = values[:, index_start : index_end + 1, ...]
    analysed_time = time_values[index_start : index_end + 1]
    if analysed_time.size < 2:
        raise ValueError("the selected time_window must contain at least two samples")

    # ---------------------------- PCA computation ----------------------------
    print("Computing PCA")

    if analysed_data.ndim == 4:
        data_for_pca = np.mean(analysed_data, axis=3)
    else:
        data_for_pca = analysed_data

    if data_for_pca.ndim == 3:
        averaged_data = np.mean(data_for_pca, axis=2)
    else:
        averaged_data = data_for_pca

    activation_patterns, variance = _compute_pca(averaged_data)

    # ----------------------- PCA on randomized data --------------------------
    variance_permutations: FloatArray | None = None
    if permutations > 0:
        print("Computing PCA on randomized data")
        generator = np.random.default_rng(
            None if random_state is None else int(random_state)
        )
        randomized_variance = np.empty((variance.size, permutations), dtype=float)

        for permutation in range(permutations):
            randomized_data = _randomize_data(
                averaged_data,
                randomization,
                generator,
            )
            _, randomized_variance[:, permutation] = _compute_pca(randomized_data)
            print(f"Permutation number {permutation + 1} / {permutations}")

        variance_permutations = np.mean(randomized_variance, axis=1)
        maximum_randomized_variance = float(np.max(randomized_variance))
        significant_indices = np.flatnonzero(variance > maximum_randomized_variance)
        significant_networks = significant_indices.astype(np.int64) + 1
        number_time_series = significant_networks.size

        if number_time_series == 0:
            print("There are no PCs (brain networks) that survived MCS")
        else:
            print("Percentage of variance explained by significant PCs (brain networks)")
            print(variance[significant_indices])
    else:
        significant_networks = "MCS has not been run"
        number_time_series = variance.size

    # ---------------------- Normalize eigenvector signs ----------------------
    activation_patterns = _normalize_eigenvector_signs(
        activation_patterns,
        sign_method,
    )

    # -------------------------- Prepare output -------------------------------
    print("Preparing output")
    if number_time_series == 0:
        time_series = np.empty(0, dtype=float)
    else:
        time_series = _network_time_series(
            analysed_data,
            activation_patterns,
            number_time_series,
        )

    return BROADNESSResult(
        Variance_BrainNetworks=np.asarray(variance, dtype=float),
        Significant_BrainNetworks=significant_networks,
        ActivationPatterns_BrainNetworks=np.asarray(
            activation_patterns,
            dtype=float,
        ),
        TimeSeries_BrainNetworks=time_series,
        Time=np.asarray(analysed_time, dtype=float),
        OriginalData=analysed_data,
        VariancePermutations=variance_permutations,
    )


__all__ = ["BROADNESSResult", "BROADNESS_NetworkEstimation"]
