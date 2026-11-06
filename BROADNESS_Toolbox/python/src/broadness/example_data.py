"""Load the published Zenodo example data for the guided Python scripts.

This module is a convenience for the three examples, not part of the analysis
algorithm. The tutorial subset is deliberately unsuitable for inference.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Sequence

import hdf5storage
import numpy as np
from numpy.typing import NDArray
from scipy.io import loadmat


AVERAGED_NAME = "DataReduced_AveragedOverParticipants_Example.mat"
DATASET_URL = "https://zenodo.org/records/17048137"


@dataclass(frozen=True)
class PublishedExampleData:
    data: NDArray[np.number]
    time: NDArray[np.float64]
    mni_coords: NDArray[np.float64]
    participant_files: tuple[Path, ...]
    source_indices: NDArray[np.int64]
    time_indices: NDArray[np.int64]
    tutorial_mode: bool


def default_data_path() -> Path:
    """Return the repository's top-level Data directory when found."""
    candidates = [Path.cwd(), *Path.cwd().parents, *Path(__file__).resolve().parents]
    for candidate in candidates:
        if (candidate / "BROADNESS_Toolbox" / "python").is_dir():
            return candidate / "Data"
    return Path.cwd() / "Data"


def _indices(length: int, limit: int | None) -> NDArray[np.int64]:
    if limit is None or limit >= length:
        return np.arange(length, dtype=np.int64)
    if limit < 2:
        raise ValueError("tutorial source/time limits must be at least 2")
    return np.unique(np.linspace(0, length - 1, limit, dtype=np.int64))


def load_published_example(
    data_path: str | Path | None = None,
    *,
    participants: bool = False,
    participant_files: Sequence[str] | None = None,
    max_participants: int | None = 6,
    tutorial_mode: bool = True,
    tutorial_sources: int = 128,
    tutorial_timepoints: int = 160,
    allow_large: bool = False,
) -> PublishedExampleData:
    """Load the averaged file or selected participant files from Zenodo.

    Files are expected directly in ``Data`` at the repository root. When
    ``participants=True``, the averaged file supplies its 776-sample time
    vector and participant data are trimmed to that common duration, as in
    the MATLAB examples. The tutorial subset samples sources and time across
    the full range; it is for learning the workflow only.
    """
    directory = default_data_path() if data_path is None else Path(data_path).expanduser()
    averaged_file = directory / AVERAGED_NAME
    if not averaged_file.is_file():
        raise FileNotFoundError(
            f"Missing {averaged_file}. Download the averaged .mat file from {DATASET_URL} "
            "and put it directly in the top-level Data folder."
        )
    if not isinstance(tutorial_mode, bool):
        raise TypeError("tutorial_mode must be True or False")
    if not isinstance(allow_large, bool):
        raise TypeError("allow_large must be True or False")
    if max_participants is not None and (not isinstance(max_participants, int) or max_participants < 1):
        raise ValueError("max_participants must be a positive integer or None")

    averaged = hdf5storage.loadmat(averaged_file, variable_names=["time", "data"] if not participants else ["time"])
    time = np.asarray(averaged["time"], dtype=float).reshape(-1)
    if time.size < 2 or not np.all(np.isfinite(time)) or np.any(np.diff(time) <= 0):
        raise ValueError("the averaged file must contain a finite, increasing time vector")

    if participants:
        if participant_files is None:
            selected_files = sorted(directory.glob("SUBJ*.mat"))
            if max_participants is not None:
                selected_files = selected_files[:max_participants]
        else:
            selected_files = [directory / name for name in participant_files]
        if not selected_files:
            raise FileNotFoundError(
                f"No SUBJ*.mat files found in {directory}. Download a participant ZIP from {DATASET_URL}, "
                "extract a few .mat files, and put them directly in this folder."
            )
        missing = [path.name for path in selected_files if not path.is_file()]
        if missing:
            raise FileNotFoundError(f"Missing participant files in {directory}: {', '.join(missing)}")
        if len({path.name for path in selected_files}) != len(selected_files):
            raise ValueError("participant_files contains duplicates")
        if len(selected_files) > 20 and not allow_large:
            raise MemoryError(
                f"Reading {len(selected_files)} participant files can take substantial time "
                "and disk I/O, even in tutorial mode. Select at most 20 files or set "
                "ALLOW_LARGE = True in the example after checking available resources."
            )
        first = loadmat(selected_files[0], variable_names=["Data"])
        if "Data" not in first:
            raise ValueError(f"{selected_files[0].name} has no Data variable")
        shape = np.asarray(first["Data"]).shape
        if len(shape) != 3 or shape[1] < time.size:
            raise ValueError("participant Data must be sources × time × conditions and cover the averaged time vector")
        n_sources, _, n_conditions = shape
    else:
        if "data" not in averaged:
            raise ValueError(f"{AVERAGED_NAME} has no data variable")
        shape = np.asarray(averaged["data"]).shape
        if len(shape) != 3 or shape[1] != time.size:
            raise ValueError("averaged data must be sources × time × conditions")
        n_sources, _, n_conditions = shape
        selected_files = []

    source_indices = _indices(n_sources, tutorial_sources if tutorial_mode else None)
    time_indices = _indices(time.size, tutorial_timepoints if tutorial_mode else None)
    estimated_bytes = (len(source_indices) * len(time_indices) * n_conditions *
                       (len(selected_files) if participants else 1) * 8)
    if estimated_bytes > 2 * 1024**3 and not allow_large:
        raise MemoryError(
            f"The selected data array needs about {estimated_bytes / 1024**3:.1f} GiB, "
            "before PCA and statistics. Select fewer participants or set "
            "ALLOW_LARGE = True in the example after checking available memory."
        )
    coordinate_file = files("broadness").joinpath("assets/MNI152_8mm_coord_dyi.mat")
    coordinates = np.asarray(loadmat(coordinate_file)["MNI8"], dtype=float)
    if coordinates.shape != (n_sources, 3):
        raise ValueError("the published data source count does not match the bundled MNI coordinates")
    coordinates = coordinates[source_indices]

    if participants:
        data = np.empty((len(source_indices), len(time_indices), n_conditions, len(selected_files)), dtype=np.float64)
        for index, path in enumerate(selected_files):
            values = np.asarray(first["Data"] if index == 0 else loadmat(path, variable_names=["Data"])["Data"])
            if values.ndim != 3 or values.shape[0] != n_sources or values.shape[2] != n_conditions or values.shape[1] < time.size:
                raise ValueError(f"{path.name} has dimensions inconsistent with the first participant")
            data[:, :, :, index] = values[np.ix_(source_indices, time_indices, np.arange(n_conditions))]
            print(f"Loaded participant {index + 1}/{len(selected_files)}: {path.name}")
            if index == 0:
                del first
    else:
        data = np.asarray(averaged["data"])[np.ix_(source_indices, time_indices, np.arange(n_conditions))]

    print("Tutorial subset (software learning only; no scientific inference)"
          if tutorial_mode else "Full spatial and temporal data")
    print(f"Data shape (sources, time, conditions{', participants' if participants else ''}): {data.shape}")
    return PublishedExampleData(
        data=data, time=time[time_indices], mni_coords=coordinates,
        participant_files=tuple(selected_files), source_indices=source_indices,
        time_indices=time_indices, tutorial_mode=tutorial_mode,
    )


__all__ = ["PublishedExampleData", "default_data_path", "load_published_example"]
