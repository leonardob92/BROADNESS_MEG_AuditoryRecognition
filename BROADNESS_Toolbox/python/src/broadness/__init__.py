"""Python implementation of the BROADNESS toolbox.

The public package preserves the modular structure and established function
names of the MATLAB BROADNESS toolbox.
"""

from .BROADNESS_Startup import BROADNESSPaths, BROADNESS_Startup
from .BROADNESS_EffectiveDimensionality import BROADNESS_EffectiveDimensionality
from .BROADNESS_NetworkEstimation import (
    BROADNESSResult,
    BROADNESS_NetworkEstimation,
)
from .BROADNESS_AlternativeNetworkEstimation_ICA import (
    BROADNESSICAResult,
    BROADNESS_AlternativeNetworkEstimation_ICA,
)
from .BROADNESS_Visualizer import BROADNESSVisualization, BROADNESS_Visualizer
from .BROADNESS_FDRCorrection import BROADNESS_FDRCorrection
from .BROADNESS_BehaviourCorrelation import BROADNESS_BehaviourCorrelation
from .BROADNESS_RepeatedMeasuresANOVA import (
    BROADNESSANOVAResult,
    BROADNESS_RepeatedMeasuresANOVA,
)
from .BROADNESS_SignificantSamplesToWindows import BROADNESS_SignificantSamplesToWindows
from ._temporal_clusters import TemporalClusterSettings, TemporalClusterStats
from .BROADNESS_TemporalClusterCorrection import BROADNESS_TemporalClusterCorrection
from .BROADNESS_TemporalClusterCorrection_ANOVA import BROADNESS_TemporalClusterCorrection_ANOVA
from .BROADNESS_TemporalClusterCorrection_Groups import BROADNESS_TemporalClusterCorrection_Groups
from .BROADNESS_TemporalClusterCorrection_GroupANOVA import BROADNESS_TemporalClusterCorrection_GroupANOVA
from .BROADNESS_PhaseSpace_RQA import BROADNESS_PhaseSpace_RQA
from .BROADNESS_PhaseSpaceStatistics import BROADNESS_PhaseSpaceStatistics
from .BROADNESS_SpatialActivationClustering import BROADNESS_SpatialActivationClustering
from .BROADNESS_SpatialPatternStatistics import BROADNESS_SpatialPatternStatistics
from .BROADNESS_SpatialPatternGroupStatistics import BROADNESS_SpatialPatternGroupStatistics

__version__ = "0.1.0"

__all__ = [
    "BROADNESSPaths",
    "BROADNESSANOVAResult",
    "BROADNESSICAResult",
    "BROADNESSResult",
    "BROADNESS_AlternativeNetworkEstimation_ICA",
    "BROADNESS_BehaviourCorrelation",
    "BROADNESS_FDRCorrection",
    "BROADNESS_EffectiveDimensionality",
    "BROADNESS_NetworkEstimation",
    "BROADNESS_PhaseSpace_RQA",
    "BROADNESS_PhaseSpaceStatistics",
    "BROADNESS_RepeatedMeasuresANOVA",
    "BROADNESS_SignificantSamplesToWindows",
    "BROADNESS_SpatialActivationClustering",
    "BROADNESS_SpatialPatternStatistics",
    "BROADNESS_SpatialPatternGroupStatistics",
    "BROADNESS_TemporalClusterCorrection",
    "BROADNESS_TemporalClusterCorrection_ANOVA",
    "BROADNESS_TemporalClusterCorrection_GroupANOVA",
    "BROADNESS_TemporalClusterCorrection_Groups",
    "BROADNESS_Startup",
    "BROADNESSVisualization",
    "BROADNESS_Visualizer",
    "TemporalClusterSettings",
    "TemporalClusterStats",
    "__version__",
]
