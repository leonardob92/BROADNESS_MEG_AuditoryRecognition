# %% 0) STARTUP AND USER SETTINGS
"""General BROADNESS example on the published Zenodo data.

Run cells in order in VS Code. Tutorial mode samples the real data for a quick
software lesson; set TUTORIAL_MODE = False only for a planned research run.
See docs/tutorials/01_GENERAL.md for a step-by-step explanation.
"""

from pathlib import Path

from broadness import (
    BROADNESS_AlternativeNetworkEstimation_ICA,
    BROADNESS_EffectiveDimensionality,
    BROADNESS_NetworkEstimation,
    BROADNESS_PhaseSpace_RQA,
    BROADNESS_PhaseSpaceStatistics,
    BROADNESS_SpatialActivationClustering,
    BROADNESS_Visualizer,
)
from broadness.example_data import default_data_path, load_published_example

data_path = default_data_path()  # Change to Path("/your/data/folder") if needed.
output_path = default_data_path().parent / "BROADNESS_Toolbox" / "Output"
TUTORIAL_MODE = False
USE_PARTICIPANTS = False  # False = 51-MB averaged file; True = participant .mat files.
MAX_PARTICIPANTS = 3       # Used only when USE_PARTICIPANTS is True.
ALLOW_LARGE = False        # Explicit opt-in if the selected array exceeds 2 GiB.
CONDITION_LABELS = ["Memorized", "NewT1", "NewT2", "NewT3", "NewT4"]
FIGURE_MODE = "show"      # "show", "save", "both", or "off".
print("Data folder:", data_path)


# %% 1) PERFORM BROADNESS (ESSENTIAL INPUTS)
example = load_published_example(
    data_path, participants=USE_PARTICIPANTS,
    max_participants=MAX_PARTICIPANTS, tutorial_mode=TUTORIAL_MODE,
    allow_large=ALLOW_LARGE,
)
data, time, mni_coords = example.data, example.time, example.mni_coords
if data.shape[2] != len(CONDITION_LABELS):
    raise ValueError("Update CONDITION_LABELS to match the number of conditions")
BROADNESS = BROADNESS_NetworkEstimation(data, time, components_num=2)
print("Network time-series shape:", BROADNESS.TimeSeries_BrainNetworks.shape)
print("Explained variance (%):", BROADNESS.Variance_BrainNetworks[:2])


# %% 1b) ALTERNATIVE PCA SETTINGS (OPTIONAL)
# Run this cell if you want to compare settings. It does not replace BROADNESS.
RUN_ALTERNATIVE_PCA = False
if RUN_ALTERNATIVE_PCA:
    BROADNESS_ALTERNATIVE = BROADNESS_NetworkEstimation(
        data, time, components_num=2, sign_eigenvect="occurrences",
        permutations_num=10, random_state=7,
    )
    print("Alternative PCA:", BROADNESS_ALTERNATIVE.Variance_BrainNetworks[:2])
else:
    print("Alternative PCA skipped. Set RUN_ALTERNATIVE_PCA = True to run it.")


# %% 1-extension) EFFECTIVE DIMENSIONALITY
effective_dimensionality = BROADNESS_EffectiveDimensionality(
    BROADNESS.Variance_BrainNetworks
)
print("Effective dimensionality:", effective_dimensionality)


# %% 2) VISUALIZE VARIANCE AND NETWORK TIME SERIES
spatial_coords = example.mni_coords
visualization = BROADNESS_Visualizer(
    BROADNESS,
    {
        "WhichPlots": [0, 0, 0, 0, 1],
        "ncomps": [1, 2],  # BROADNESS selections start at 1, like MATLAB.
        "ncomps_var": 2,
        "MNI_coords": spatial_coords,
        "Labels": CONDITION_LABELS,
        "OutputPath": output_path,
        "FigureMode": FIGURE_MODE,
        "FigureLayout": "individual",
    },
)
print("Saved figure files:", visualization.Figures.Files)

# %% 3) PHASE SPACE AND RECURRENCE QUANTIFICATION ANALYSIS
# Recurrence matrices grow quadratically with the number of time samples.
RUN_RQA = True
if RUN_RQA:
    RQA_BROADNESS, RQA_FIGURES = BROADNESS_PhaseSpace_RQA(
        BROADNESS, principalcomps=(1, 2), threshold=0.1,
        figuremode="off", video="off",
    )
    print("RQA measures for participant 1:")
    print(RQA_BROADNESS.RQA_metrics[0])
else:
    print("RQA skipped.")


# %% 3b) PHASE-SPACE METRICS AND STATISTICS
if RUN_RQA:
    PHASE_STATS = BROADNESS_PhaseSpaceStatistics(
        RQA_BROADNESS, conditionnames=CONDITION_LABELS,
        statistics="on" if USE_PARTICIPANTS else "off",
        testtype="parametric",
    )
    print("Speed array (time, conditions, participants):", PHASE_STATS.Metrics.Speed.shape)
else:
    print("Run section 3 first to calculate phase-space measures.")


# %% 4) SPATIAL ACTIVATION PATTERN CLUSTERING
RUN_SPATIAL_CLUSTERING = True
if RUN_SPATIAL_CLUSTERING:
    SPATIAL_CLUSTERING_BROADNESS = BROADNESS_SpatialActivationClustering(
        BROADNESS, principalcomps=(1, 2), nclusters=(2, 3, 4),
        evalclusters=3, mni_coords=mni_coords,
        figuremode="show", OutputPath=None if TUTORIAL_MODE else output_path,
    )
    print("Spatial clustering complete. Inspect SPATIAL_CLUSTERING_BROADNESS in Variables.")
else:
    print("Spatial clustering skipped.")


# %% 5) ALTERNATIVE NETWORK ESTIMATION WITH ICA
RUN_ICA = True
if RUN_ICA:
    BROADNESS_ICA = BROADNESS_AlternativeNetworkEstimation_ICA(
        data, time, icacomps=2,
    )
    print("ICA time-series shape:", BROADNESS_ICA.TimeSeries_BrainNetworks.shape)
else:
    print("ICA skipped.")


# %% 5b) VISUALIZE ICA RESULTS
if RUN_ICA:
    ICA_VISUALIZATION = BROADNESS_Visualizer(
        BROADNESS_ICA,
        {
            "WhichPlots": [0, 0, 1, 0, 0],
            "ncomps": [1, 2],
            "Labels": CONDITION_LABELS,
            "OutputPath": output_path,
            "FigureMode": FIGURE_MODE,
        },
    )
    print("ICA figure files:", ICA_VISUALIZATION.Figures.Files)
else:
    print("Run section 5 first to visualize ICA.")
