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
| 0 | Keep `USE_PARTICIPANTS = False` and confirm the printed `Data` folder. The example now uses the full averaged file by default. | `data_path` |
| 1 | Load the averaged file and estimate two PCA networks. | `data` has shape `(3559, 776, 5)`; `BROADNESS.TimeSeries_BrainNetworks` has shape `(776, 2, 5)` |
| 1b | Optional comparison of PCA settings. Change `RUN_ALTERNATIVE_PCA` to `True` only when wanted. | `BROADNESS_ALTERNATIVE` |
| 1-extension | Calculate effective dimensionality from the PCA variance spectrum. | `effective_dimensionality` |
| 2 | Export NIfTI images and activation tables for the two networks. | `visualization.nifti_paths`, `visualization.activation_table_paths` |
| 3 | Calculate phase trajectories and RQA measures. | `RQA_BROADNESS.RQA_metrics[0]` |
| 3b | Calculate descriptive phase-space measures. Averaged data do not support participant-level inference. | `PHASE_STATS.Metrics.Speed` |
| 4 | Cluster source activation patterns. Each cluster figure shows an oblique brain view and a view from above, with a compact PC activation profile below. | `SPATIAL_CLUSTERING_BROADNESS` |
| 5 and 5b | Estimate two ICA networks and visualize their time series. | `BROADNESS_ICA`, `ICA_VISUALIZATION` |

Section 2 writes the NIfTI images and activation tables under
`BROADNESS_Toolbox/Output/`; that folder stays on your computer and is ignored
by Git. Section 4 shows the clustering figures.

## Inspect the brain from two angles

Run sections 0 and 1, then click **Run Cell** for section 4. Each cluster
figure shows an oblique view and a view from above, with the PC activation
profile below them. These are ordinary figures that fit in the Interactive
Window; there is no dragging or extra 3D renderer to install. Restart the
Interactive Window kernel once if you previously ran the old widget version.

The default example uses the full averaged array, about 55 MB in memory, and
does not read participant files. Saved PNG/PDF figures use the same two views.

## Try a few individual participants

The participant files are inside large ZIP archives on Zenodo. Extract a few
`SUBJ_*.mat` files into the same `Data/` folder **and keep the averaged file**,
which supplies the time vector. In section 0 set `USE_PARTICIPANTS = True` and
keep `MAX_PARTICIPANTS = 3`. Restart the Interactive Window and run from section
0. The script prints exactly which files it selected, in participant order.
Participant data add a fourth axis: sources × time × conditions × participants.

## When moving beyond the lesson

The example now starts with `TUTORIAL_MODE = False` and the full averaged
file. Setting it to `True` selects 128 of 3,559 sources and 160 of 776 time
samples for a quicker software lesson, but the spatial maps will be sparse
and unsuitable for scientific interpretation. When using participants, set
`MAX_PARTICIPANTS = None` only after deciding how many files and how much
memory to use. Full PCA,
participant projection, figures, and RQA can take substantially longer; RQA
also holds time × time matrices. The MATLAB and Python implementations have
not yet been directly compared numerically on the complete dataset. The loader
stops before allocating a data array larger than 2 GiB unless you explicitly
set `ALLOW_LARGE = True` after checking your computer's memory. It also stops
before reading more than 20 participant files by default, even in tutorial
mode, because those files require substantial disk I/O.
