"""ICA-based BROADNESS network estimation using the JADE algorithm.

The public function mirrors the MATLAB toolbox's input axes and output fields.
Its ICA decomposition follows the real-valued JADE procedure used by
``jadeR.m``: PCA whitening, fourth-order cumulant matrices, joint Jacobi
diagonalization, energy ordering, and first-source sign normalization.
Independent components remain identifiable only up to sign and permutation.
"""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import svd

from .BROADNESS_NetworkEstimation import _compute_pca, _numeric_data, _time_vector


FloatArray = NDArray[np.float64]


@dataclass(slots=True)
class BROADNESSICAResult:
    """MATLAB-compatible outputs of :func:`BROADNESS_AlternativeNetworkEstimation_ICA`."""

    TimeSeries_BrainNetworks: FloatArray
    ActivationPatterns_BrainNetworks: FloatArray
    numICAcomps: int
    Time: FloatArray
    OriginalData: NDArray[np.number]


def _jade_unmixing(data: FloatArray, n_components: int) -> FloatArray:
    """Return a JADE separating matrix for sources × samples data.

    The cumulant stack uses O(m^4) memory for m requested components, as in
    MATLAB's ``jadeR``. Whitening uses a thin SVD to avoid a sources × sources
    covariance eigendecomposition when there are many brain sources.
    """
    n_sources, n_samples = data.shape
    centered = data - np.mean(data, axis=1, keepdims=True)
    left_vectors, singular_values, _ = svd(
        centered, full_matrices=False, check_finite=False, lapack_driver="gesdd"
    )
    tolerance = np.finfo(float).eps * max(data.shape) * singular_values[0]
    rank = int(np.count_nonzero(singular_values > tolerance))
    if n_components > rank:
        raise ValueError(
            f"ICA needs {n_components} independent components, but the "
            f"averaged data have rank {rank}"
        )

    scales = singular_values[:n_components] / np.sqrt(n_samples)
    whitener = (left_vectors[:, :n_components] / scales).T
    inverse_whitener = left_vectors[:, :n_components] * scales
    whitened = whitener @ centered

    m = n_components
    n_cumulants = m * (m + 1) // 2
    identity = np.eye(m)
    cumulants = np.empty((n_cumulants, m, m), dtype=float)
    index = 0
    for i in range(m):
        basis_i = identity[:, i]
        cumulants[index] = (
            (whitened * whitened[i] ** 2) @ whitened.T / n_samples
            - identity
            - 2.0 * np.outer(basis_i, basis_i)
        )
        index += 1
        for j in range(i):
            basis_j = identity[:, j]
            cumulants[index] = np.sqrt(2.0) * (
                (whitened * (whitened[i] * whitened[j]))
                @ whitened.T
                / n_samples
                - np.outer(basis_i, basis_j)
                - np.outer(basis_j, basis_i)
            )
            index += 1

    _, rotation = np.linalg.eigh(cumulants[0])
    cumulants = rotation.T @ cumulants @ rotation

    threshold = 1.0 / (100.0 * np.sqrt(n_samples))
    for _ in range(100):
        rotated = False
        for p in range(m - 1):
            for q in range(p + 1, m):
                g1 = cumulants[:, p, p] - cumulants[:, q, q]
                g2 = cumulants[:, p, q] + cumulants[:, q, p]
                g11 = float(g1 @ g1)
                g22 = float(g2 @ g2)
                g12 = float(g1 @ g2)
                on = g11 - g22
                off = 2.0 * g12
                angle = 0.5 * np.arctan2(off, on + np.hypot(on, off))
                if abs(angle) <= threshold:
                    continue

                rotated = True
                cosine, sine = np.cos(angle), np.sin(angle)
                givens = np.array([[cosine, -sine], [sine, cosine]])
                rotation[:, [p, q]] = rotation[:, [p, q]] @ givens
                cumulants[:, [p, q], :] = (
                    givens.T @ cumulants[:, [p, q], :]
                )
                cumulants[:, :, [p, q]] = cumulants[:, :, [p, q]] @ givens
        if not rotated:
            break
    else:
        warnings.warn(
            "JADE did not converge within 100 sweeps; inspect the ICA result",
            RuntimeWarning,
            stacklevel=2,
        )

    separating = rotation.T @ whitener
    mixing = inverse_whitener @ rotation
    order = np.argsort(np.sum(mixing * mixing, axis=0))[::-1]
    separating = separating[order]
    separating *= np.where(separating[:, :1] < 0.0, -1.0, 1.0)
    assert separating.shape == (m, n_sources)
    return np.asarray(separating, dtype=float)


def _project_ica(data: NDArray[np.number], separating: FloatArray) -> FloatArray:
    """Apply the common JADE separating matrix to each original data slice."""
    if data.ndim == 2:
        return np.asarray((separating @ data).T, dtype=float)
    if data.ndim == 3:
        return np.asarray(np.einsum("ks,stc->tkc", separating, data), dtype=float)
    return np.asarray(np.einsum("ks,stcp->tkcp", separating, data), dtype=float)


def BROADNESS_AlternativeNetworkEstimation_ICA(
    data: ArrayLike,
    time: ArrayLike,
    *,
    icacomps: int | None = None,
    total_varexp: float = 95.0,
) -> BROADNESSICAResult:
    """Estimate independent brain networks with MATLAB-compatible JADE.

    Data have shape sources × time × [conditions] × [participants]. A common
    decomposition is fitted to the mean across conditions and participants.
    ``icacomps`` selects the number of components directly; otherwise the
    smallest number explaining at least ``total_varexp`` percent of PCA
    variance is used. The unmixing matrix is applied to each original slice.
    """
    values = _numeric_data(data)
    time_values = _time_vector(time, values.shape[1])
    if (
        isinstance(total_varexp, (bool, np.bool_))
        or not isinstance(total_varexp, (int, float, np.integer, np.floating))
    ):
        raise TypeError("total_varexp must be a real number in (0, 100]")
    try:
        variance_threshold = float(total_varexp)
    except (TypeError, ValueError) as error:
        raise TypeError("total_varexp must be a real number in (0, 100]") from error
    if not np.isfinite(variance_threshold) or not 0 < variance_threshold <= 100:
        raise ValueError("total_varexp must be in (0, 100]")

    if values.ndim == 4:
        averaged_data = np.mean(values, axis=(2, 3), dtype=float)
    elif values.ndim == 3:
        averaged_data = np.mean(values, axis=2, dtype=float)
    else:
        averaged_data = np.asarray(values, dtype=float)

    max_components = min(averaged_data.shape)
    if icacomps is not None:
        if isinstance(icacomps, (bool, np.bool_)) or not isinstance(
            icacomps, (int, np.integer)
        ):
            raise TypeError("icacomps must be a positive integer")
        n_components = int(icacomps)
        if not 1 <= n_components <= max_components:
            raise ValueError(
                f"icacomps must be between 1 and {max_components}"
            )
    else:
        _, explained_variance = _compute_pca(averaged_data)
        cumulative = np.cumsum(explained_variance)
        n_components = int(np.searchsorted(cumulative, variance_threshold) + 1)
        n_components = min(n_components, explained_variance.size)

    separating = _jade_unmixing(averaged_data, n_components)
    time_series = _project_ica(values, separating)
    covariance = np.cov(averaged_data)
    activation_patterns = (separating @ covariance).T
    return BROADNESSICAResult(
        TimeSeries_BrainNetworks=time_series,
        ActivationPatterns_BrainNetworks=np.asarray(activation_patterns, dtype=float),
        numICAcomps=n_components,
        Time=time_values,
        OriginalData=values,
    )


__all__ = ["BROADNESSICAResult", "BROADNESS_AlternativeNetworkEstimation_ICA"]
