# Tutorial 1: estimate networks from the published data

This follows the section order of MATLAB's `BROADNESS_Example.m`. The script is
[`01_BROADNESS_Example.py`](../../examples/01_BROADNESS_Example.py). It uses
the [published Zenodo dataset](https://zenodo.org/records/17048137).

## Before opening the script

1. Complete the [VS Code first steps](../VS_CODE_FIRST_STEPS.md), including
   selecting this project's `.venv` interpreter.
2. Download `DataReduced_AveragedOverParticipants_Example.mat` from Zenodo.
   Place it directly in `BROADNESS_MEG_AuditoryRecognition/Data/`, beside
   `BROADNESS_Toolbox/`. See the [data-folder diagram](../VS_CODE_FIRST_STEPS.md#6-download-real-example-data-optional).
3. Open the script in VS Code. Use **Run Cell** above each `# %%` section,
   starting at section 0. Keep the same Interactive Window open.

## Walk through the sections

| Section | What to do | What to inspect in Variables |
| --- | --- | --- |
| 0 | Keep `TUTORIAL_MODE = True` and `USE_PARTICIPANTS = False`. Confirm the printed `Data` folder. | `data_path` |
| 1 | Load the averaged file and estimate two PCA networks. | `data` has shape `(128, 160, 5)`; `BROADNESS.TimeSeries_BrainNetworks` has shape `(160, 2, 5)` |
| 1b | Optional comparison of PCA settings. Change `RUN_ALTERNATIVE_PCA` to `True` only when wanted. | `BROADNESS_ALTERNATIVE` |
| 1-extension | Calculate effective dimensionality from the PCA variance spectrum. | `effective_dimensionality` |
| 2 | Display variance and network time-series figures. | `visualization` and its `Figures` |
| 2b | Click **Run Cell** to plot the 3D map. This cell uses the full averaged file, so the map has all 3,559 source locations. | `BROADNESS_SPATIAL`, `spatial_visualization` |
| 3 | Calculate phase trajectories and RQA measures. | `RQA_BROADNESS.RQA_metrics[0]` |
| 3b | Calculate descriptive phase-space measures. Averaged data do not support participant-level inference. | `PHASE_STATS.Metrics.Speed` |
| 4 | Cluster source activation patterns. | `SPATIAL_CLUSTERING_BROADNESS` |
| 5 and 5b | Estimate two ICA networks and visualize their time series. | `BROADNESS_ICA`, `ICA_VISUALIZATION` |

In section 0, change `FIGURE_MODE` to `"save"` to write figures under
`BROADNESS_Toolbox/Output/`. That folder stays on your computer and is ignored
by Git. `"off"` suppresses figures.

## Rotate the 3D brain in VS Code

1. Run sections 0 and 1, then click **Run Cell** for section 2b. Keep
   `FIGURE_MODE = "show"` (or `"both"`). The visualizer enables the interactive
   widget automatically in the VS Code Interactive Window.
2. Drag inside the resulting 3D figure to rotate it; use the mouse wheel or
   trackpad scroll to zoom. The ordinary inline image and saved PNG are static.

If the figure stays static, run `bash setup_vscode.sh` from the `python/`
folder on macOS/Linux, or `.\setup_vscode.bat` on Windows. Then restart the
**Interactive Window kernel** and rerun sections 0, 1, and 2b. The setup script
installs the small `ipympl` plotting widget into this project's `.venv`.

Section 2b uses the complete averaged file by default, even though sections 1
and 2 use the small tutorial subset. It does this so the spatial pattern is
based on the same spatial and temporal resolution as the MATLAB example. The
full averaged array is about 55 MB in memory; it does not read participant
files. `BROADNESS_SPATIAL` is separate from the earlier `BROADNESS` result.
The 3D widget may rotate less smoothly than MATLAB because Matplotlib redraws
the brain surface as you drag.

## Try a few individual participants

The participant files are inside large ZIP archives on Zenodo. Extract a few
`SUBJ_*.mat` files into the same `Data/` folder **and keep the averaged file**,
which supplies the time vector. In section 0 set `USE_PARTICIPANTS = True` and
keep `MAX_PARTICIPANTS = 3`. Restart the Interactive Window and run from section
0. The script prints exactly which files it selected, in participant order.
Participant data add a fourth axis: sources × time × conditions × participants.

## When moving beyond the lesson

`TUTORIAL_MODE = True` deliberately selects 128 of 3,559 sources and 160 of
776 time samples across the real dataset. It is a software lesson, not a
reproduction of the paper's analysis. For full resolution, set
`TUTORIAL_MODE = False`; when using participants, set `MAX_PARTICIPANTS = None`
only after deciding how many files and how much memory to use. Full PCA,
participant projection, figures, and RQA can take substantially longer; RQA
also holds time × time matrices. The MATLAB and Python implementations have
not yet been directly compared numerically on the complete dataset. The loader
stops before allocating a data array larger than 2 GiB unless you explicitly
set `ALLOW_LARGE = True` after checking your computer's memory. It also stops
before reading more than 20 participant files by default, even in tutorial
mode, because those files require substantial disk I/O.
