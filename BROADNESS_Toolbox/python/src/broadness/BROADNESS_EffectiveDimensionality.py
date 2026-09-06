"""Effective dimensionality estimation for BROADNESS.

=========================================================================
 USER DOCUMENTATION
=========================================================================

This module provides the public function
``BROADNESS_EffectiveDimensionality``. It estimates how many principal
components (brain networks) contribute meaningfully to an eigenspectrum.

INPUT
-----
``eigenspectrum`` contains the variance explained by the components. For a
standard BROADNESS analysis, use::

    BROADNESS.Variance_BrainNetworks

The first dimension must contain the components. A one-dimensional vector
returns one effective-dimensionality value. Additional dimensions are
supported, and one value is computed independently for each eigenspectrum.

EXAMPLE
-------
Run the analysis using::

    from broadness import BROADNESS_EffectiveDimensionality

    ED = BROADNESS_EffectiveDimensionality(
        BROADNESS.Variance_BrainNetworks
    )

OUTPUT
------
``ED`` is the estimated number of components to retain. It is computed as
the rounded participation ratio::

    (sum(eigenspectrum) ** 2) / sum(eigenspectrum ** 2)

An eigenspectrum concentrated in one component therefore approaches an
effective dimensionality of one, whereas an eigenspectrum distributed
equally across several components approaches the number of those components.

If you use this toolbox, please cite the first BROADNESS paper:
Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
F., Testa, C., Vuust, P., Kringelbach, M. L., & Rosso, M. (2025).
BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in
Auditory Memory Networks. Advanced Science.
https://doi.org/10.1002/advs.202507878

Toolbox authors
---------------
Mattia Rosso - Center for Music in the Brain, Aarhus University -
mattia.rosso@clin.au.dk
Chiara Malvaso - Department of Physics, University of Bologna -
chiara.malvaso@studio.unibo.it
Leonardo Bonetti - Center for Music in the Brain, Aarhus University; Centre
for Eudaimonia and Human Flourishing, Linacre College, University of Oxford -
leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


# =========================================================================
#  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS)
#  EFFECTIVE DIMENSIONALITY ESTIMATION
# =========================================================================

IntegerArray = NDArray[np.int64]


def BROADNESS_EffectiveDimensionality(
    eigenspectrum: ArrayLike,
) -> int | IntegerArray:
    """Estimate the effective dimensionality of one or more eigenspectra.

    Parameters
    ----------
    eigenspectrum
        Variance explained by the components. Components occupy the first
        dimension. A one-dimensional input is treated as one eigenspectrum.

    Returns
    -------
    int or numpy.ndarray
        Rounded effective dimensionality. A single eigenspectrum returns an
        integer; additional input dimensions return an integer array after
        dimensions of length one are removed, matching MATLAB ``squeeze``.

    Notes
    -----
    Effective dimensionality is the participation ratio
    ``(sum(lambda) ** 2) / sum(lambda ** 2)``. The result is invariant to a
    common positive scaling of the eigenspectrum.
    """
    print("Computing effective dimensionality")

    values = np.asarray(eigenspectrum)

    if values.ndim == 0:
        raise ValueError("eigenspectrum must contain at least one component")
    if values.size == 0 or values.shape[0] == 0:
        raise ValueError("eigenspectrum cannot be empty")
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("eigenspectrum must contain real numeric values")

    values = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError("eigenspectrum cannot contain NaN or infinite values")
    if np.any(values < 0):
        raise ValueError("eigenspectrum cannot contain negative values")

    numerator = np.sum(values, axis=0) ** 2
    denominator = np.sum(values**2, axis=0)

    # Match the realmin safeguard in the MATLAB implementation.
    denominator = np.maximum(denominator, np.finfo(np.float64).tiny)

    # MATLAB rounds positive half-integers away from zero. NumPy's round uses
    # ties-to-even, so floor(x + 0.5) is used for MATLAB-compatible results.
    effective_dimensionality = np.floor(
        numerator / denominator + 0.5
    ).astype(np.int64)
    effective_dimensionality = np.squeeze(effective_dimensionality)

    if effective_dimensionality.ndim == 0:
        return int(effective_dimensionality)
    return effective_dimensionality


__all__ = ["BROADNESS_EffectiveDimensionality"]
