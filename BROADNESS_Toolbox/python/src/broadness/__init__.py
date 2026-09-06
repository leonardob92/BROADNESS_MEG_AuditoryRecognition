"""Python implementation of the BROADNESS toolbox.

The public package preserves the modular structure and established function
names of the MATLAB BROADNESS toolbox.
"""

from .BROADNESS_Startup import BROADNESSPaths, BROADNESS_Startup

__version__ = "0.1.0"

__all__ = ["BROADNESSPaths", "BROADNESS_Startup", "__version__"]

