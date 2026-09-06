"""Initialize the Python BROADNESS toolbox and locate shared resources.

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
from pathlib import Path
import warnings


# =========================================================================
#  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS)
#  STARTUP
# =========================================================================


@dataclass(frozen=True, slots=True)
class BROADNESSPaths:
    """Resolved BROADNESS directories and external resource locations."""

    path_home: Path
    functions: Path
    external: Path
    nifti_tools: Path
    mni_coordinates: Path | None
    cerebellum_coordinates: Path | None
    brain_mask: Path | None
    nifti_template: Path | None
    matlab_brain_template: Path | None


def _toolbox_candidate(path: Path) -> Path | None:
    """Return the toolbox root represented by ``path``, when present."""
    if (path / "BROADNESS_Functions").is_dir() and (
        path / "BROADNESS_External"
    ).is_dir():
        return path

    nested_toolbox = path / "BROADNESS_Toolbox"
    if (nested_toolbox / "BROADNESS_Functions").is_dir() and (
        nested_toolbox / "BROADNESS_External"
    ).is_dir():
        return nested_toolbox

    return None


def _resolve_toolbox_root(path_home: str | Path | None) -> Path:
    """Resolve an explicitly supplied or automatically discovered toolbox."""
    if path_home is not None:
        if not isinstance(path_home, (str, Path)):
            raise TypeError("path_home must be a string or pathlib.Path")

        requested_path = Path(path_home).expanduser().resolve()
        if not requested_path.is_dir():
            raise FileNotFoundError(
                f"BROADNESS toolbox directory not found: {requested_path}"
            )

        root = _toolbox_candidate(requested_path)
        if root is None:
            raise FileNotFoundError(
                "The selected directory does not contain the expected "
                "BROADNESS_Functions and BROADNESS_External folders: "
                f"{requested_path}"
            )
        return root

    candidates = [Path.cwd(), *Path(__file__).resolve().parents]
    for candidate in candidates:
        root = _toolbox_candidate(candidate)
        if root is not None:
            return root.resolve()

    raise FileNotFoundError(
        "Could not locate the BROADNESS toolbox. Pass the BROADNESS_Toolbox "
        "directory to BROADNESS_Startup(path_home)."
    )


def _resource_path(external: Path, filename: str) -> Path | None:
    """Return an optional shared resource and warn when it is unavailable."""
    resource = external / filename
    if resource.is_file():
        return resource

    warnings.warn(
        f"Optional BROADNESS resource not found: {resource}",
        UserWarning,
        stacklevel=3,
    )
    return None


def BROADNESS_Startup(path_home: str | Path | None = None) -> BROADNESSPaths:
    """Initialize BROADNESS and locate its shared external resources.

    Parameters
    ----------
    path_home
        Path to ``BROADNESS_Toolbox``. The repository root containing that
        folder is also accepted. When omitted, the function searches the
        current directory and the installed package's parent directories.

    Returns
    -------
    BROADNESSPaths
        Resolved toolbox directories and paths to the available shared MNI,
        NIfTI, and MATLAB visualization resources.

    Notes
    -----
    Unlike MATLAB startup, this function does not modify global import paths.
    Installing the Python package provides that functionality without runtime
    side effects. Participant data are never searched for or loaded here.
    """
    print("Checking BROADNESS toolbox directories")

    root = _resolve_toolbox_root(path_home)
    functions = root / "BROADNESS_Functions"
    external = root / "BROADNESS_External"
    nifti_tools = external / "NIfTI_20140122"

    if not nifti_tools.is_dir():
        raise FileNotFoundError(
            f"BROADNESS NIfTI tools directory not found: {nifti_tools}"
        )

    paths = BROADNESSPaths(
        path_home=root,
        functions=functions,
        external=external,
        nifti_tools=nifti_tools,
        mni_coordinates=_resource_path(external, "MNI152_8mm_coord_dyi.mat"),
        cerebellum_coordinates=_resource_path(external, "cerebellum_coords.mat"),
        brain_mask=_resource_path(external, "MNI152_8mm_brain_diy.nii.gz"),
        nifti_template=_resource_path(external, "MNI152_T1_8mm_Template.nii.gz"),
        matlab_brain_template=_resource_path(external, "BrainTemplate_GT.fig"),
    )

    print("BROADNESS successfully initialized")
    print(f"Base directory: {root}")
    return paths


__all__ = ["BROADNESSPaths", "BROADNESS_Startup"]

