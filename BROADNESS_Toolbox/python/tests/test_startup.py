from pathlib import Path

import pytest

from broadness import BROADNESSPaths, BROADNESS_Startup


def _make_minimal_toolbox(root: Path) -> Path:
    toolbox = root / "BROADNESS_Toolbox"
    (toolbox / "BROADNESS_Functions").mkdir(parents=True)
    (toolbox / "BROADNESS_External" / "NIfTI_20140122").mkdir(parents=True)
    return toolbox


def test_startup_locates_the_complete_toolbox(capsys):
    toolbox = Path(__file__).resolve().parents[2]

    paths = BROADNESS_Startup(toolbox)

    output = capsys.readouterr().out
    assert "Checking BROADNESS toolbox directories" in output
    assert "BROADNESS successfully initialized" in output
    assert f"Base directory: {toolbox}" in output

    assert isinstance(paths, BROADNESSPaths)
    assert paths.path_home == toolbox
    assert paths.functions == toolbox / "BROADNESS_Functions"
    assert paths.external == toolbox / "BROADNESS_External"
    assert paths.nifti_tools == paths.external / "NIfTI_20140122"
    assert paths.mni_coordinates is not None
    assert paths.cerebellum_coordinates is not None
    assert paths.brain_mask is not None
    assert paths.nifti_template is not None
    assert paths.matlab_brain_template is not None


def test_startup_accepts_the_repository_root():
    toolbox = Path(__file__).resolve().parents[2]
    repository = toolbox.parent

    paths = BROADNESS_Startup(repository)

    assert paths.path_home == toolbox


def test_startup_can_discover_the_toolbox(monkeypatch):
    toolbox = Path(__file__).resolve().parents[2]
    monkeypatch.chdir(toolbox)

    paths = BROADNESS_Startup()

    assert paths.path_home == toolbox


def test_startup_rejects_an_invalid_directory(tmp_path):
    with pytest.raises(FileNotFoundError, match="does not contain"):
        BROADNESS_Startup(tmp_path)


def test_startup_rejects_an_invalid_path_type():
    with pytest.raises(TypeError, match="string or pathlib.Path"):
        BROADNESS_Startup(5)


def test_startup_requires_the_nifti_tools_directory(tmp_path):
    toolbox = tmp_path / "BROADNESS_Toolbox"
    (toolbox / "BROADNESS_Functions").mkdir(parents=True)
    (toolbox / "BROADNESS_External").mkdir()

    with pytest.raises(FileNotFoundError, match="NIfTI tools directory"):
        BROADNESS_Startup(toolbox)


def test_startup_warns_about_missing_optional_resources(tmp_path):
    toolbox = _make_minimal_toolbox(tmp_path)

    with pytest.warns(UserWarning) as warnings_record:
        paths = BROADNESS_Startup(toolbox)

    assert len(warnings_record) == 5
    assert paths.mni_coordinates is None
    assert paths.cerebellum_coordinates is None
    assert paths.brain_mask is None
    assert paths.nifti_template is None
    assert paths.matlab_brain_template is None

