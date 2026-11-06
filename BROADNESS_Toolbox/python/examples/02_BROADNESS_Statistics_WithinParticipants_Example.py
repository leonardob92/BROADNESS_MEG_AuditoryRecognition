# %% 0) STARTUP AND USER SETTINGS
"""Paired condition statistics on published participant data.

Run cells in order. Tutorial mode samples real participant data and uses FDR;
its results are for learning, not scientific inference. Read
docs/tutorials/02_WITHIN_PARTICIPANTS.md before a full-data analysis.
"""

from itertools import combinations
from pathlib import Path
from types import SimpleNamespace as Result

import numpy as np
from scipy.stats import ttest_rel

from broadness import (
    BROADNESS_BehaviourCorrelation,
    BROADNESS_FDRCorrection,
    BROADNESS_NetworkEstimation,
    BROADNESS_PhaseSpace_RQA,
    BROADNESS_PhaseSpaceStatistics,
    BROADNESS_Plot_ActivationTimeseries,
    BROADNESS_RepeatedMeasuresANOVA,
    BROADNESS_SignificantSamplesToWindows,
    BROADNESS_SpatialPatternStatistics,
    BROADNESS_TemporalClusterCorrection,
    BROADNESS_TemporalClusterCorrection_ANOVA,
)
from broadness.example_data import default_data_path, load_published_example

data_path = default_data_path()  # Change to Path("/your/data/folder") if needed.
output_path = default_data_path().parent / "BROADNESS_Toolbox" / "Output" / "Stats"
TUTORIAL_MODE = True
MAX_PARTICIPANTS = 4  # First four SUBJ*.mat files; inspect the printed file order.
ALLOW_LARGE = False   # Explicit opt-in if the selected array exceeds 2 GiB.
PARTICIPANT_FILES = None  # Or ["SUBJ_0002.mat", "SUBJ_0007.mat", ...].
SELECTED_NETWORKS = (1, 2)
CONDITION_LABELS = ["Memorized", "NewT1", "NewT2", "NewT3", "NewT4"]
CONDITION_PAIRS = ((1, 2),)  # Set to None for every condition pair.
CORRECTION = "FDR"        # "FDR" or "cluster".
ALPHA = 0.05
CLUSTER_PERMUTATIONS = 100 if TUTORIAL_MODE else 5000
FIGURE_MODE = "show"       # "show", "save", "both", or "off".


# %% 1) LOAD PARTICIPANTS AND PERFORM BROADNESS
example = load_published_example(
    data_path, participants=True, participant_files=PARTICIPANT_FILES,
    max_participants=MAX_PARTICIPANTS, tutorial_mode=TUTORIAL_MODE,
    allow_large=ALLOW_LARGE,
)
data, time, mni_coords = example.data, example.time, example.mni_coords
participant_names = [path.name for path in example.participant_files]
if data.shape[3] < 2:
    raise ValueError("This example needs at least two participant files")
if data.shape[2] != len(CONDITION_LABELS):
    raise ValueError("Update CONDITION_LABELS to match the loaded conditions")
BROADNESS = BROADNESS_NetworkEstimation(
    data, time, components_num=max(SELECTED_NETWORKS)
)
network_time_series = BROADNESS.TimeSeries_BrainNetworks
print("Participant order:", participant_names)
print("Network time-series shape:", network_time_series.shape)


# %% 2) PAIRED CONDITION TESTS ACROSS TIME
condition_pairs = (tuple(combinations(range(1, data.shape[2] + 1), 2))
                   if CONDITION_PAIRS is None else tuple(CONDITION_PAIRS))
if not condition_pairs:
    raise ValueError("Select at least one condition pair")
for first, second in condition_pairs:
    if not (1 <= first <= data.shape[2] and 1 <= second <= data.shape[2] and first < second):
        raise ValueError("CONDITION_PAIRS must use valid, increasing one-based condition numbers")
selected_series = network_time_series[:, np.array(SELECTED_NETWORKS) - 1, :, :]
n_time, n_networks, _, n_participants = selected_series.shape
p_values = np.empty((n_time, n_networks, len(condition_pairs)))
t_statistics = np.empty_like(p_values)
cohen_dz = np.empty_like(p_values)
for pair_index, (first, second) in enumerate(condition_pairs):
    difference = selected_series[:, :, first - 1, :] - selected_series[:, :, second - 1, :]
    test = ttest_rel(selected_series[:, :, first - 1, :],
                     selected_series[:, :, second - 1, :], axis=-1)
    p_values[:, :, pair_index] = test.pvalue
    t_statistics[:, :, pair_index] = test.statistic
    with np.errstate(divide="ignore", invalid="ignore"):
        cohen_dz[:, :, pair_index] = difference.mean(axis=-1) / difference.std(axis=-1, ddof=1)
cohen_dz[~np.isfinite(cohen_dz)] = np.nan

if CORRECTION.lower() == "fdr":
    significant, adjusted_p, critical_p = BROADNESS_FDRCorrection(p_values, ALPHA)
    cluster_stats = None
    condition_anova = BROADNESS_RepeatedMeasuresANOVA(selected_series, ALPHA)
elif CORRECTION.lower() == "cluster":
    significant, adjusted_p, critical_p, cluster_stats = BROADNESS_TemporalClusterCorrection(
        selected_series, condition_pairs, ALPHA,
        permutations=CLUSTER_PERMUTATIONS, clusteralpha=0.05, time=time,
    )
    anova_significant, anova_p, _, anova_details = BROADNESS_TemporalClusterCorrection_ANOVA(
        selected_series, ALPHA, permutations=CLUSTER_PERMUTATIONS,
        clusteralpha=0.05, time=time,
    )
    condition_anova = Result(Enabled=True, Significant=anova_significant,
                             ClusterCorrectedP=anova_p, ClusterResults=anova_details)
else:
    raise ValueError("CORRECTION must be 'FDR' or 'cluster'")

contrast_labels = [f"{CONDITION_LABELS[a - 1]} vs {CONDITION_LABELS[b - 1]}"
                   for a, b in condition_pairs]
CONDITION_STATS = Result(
    Settings=Result(SelectedNetworks=SELECTED_NETWORKS, ConditionLabels=CONDITION_LABELS,
                    Time=time, Alpha=ALPHA, CorrectionMethod=CORRECTION),
    PairedTests=Result(ConditionPairs=condition_pairs, ContrastLabels=contrast_labels,
                       PValues=p_values, AdjustedPValues=adjusted_p,
                       Significant=significant, TStatistics=t_statistics,
                       CohenDz=cohen_dz, ClusterResults=cluster_stats,
                       FDRCriticalP=critical_p),
    ANOVA=condition_anova,
)
print("Paired tests complete. Inspect CONDITION_STATS.PairedTests in Variables.")


# %% 3) PLOT SIGNIFICANT TIME WINDOWS
statistical_figures = []
for network_index, network in enumerate(SELECTED_NETWORKS):
    anova_windows = BROADNESS_SignificantSamplesToWindows(
        CONDITION_STATS.ANOVA.Significant[:, network_index], time
    )
    pair_windows = [BROADNESS_SignificantSamplesToWindows(
        significant[:, network_index, pair_index], time
    ) for pair_index in range(len(condition_pairs))]
    figure = BROADNESS_Plot_ActivationTimeseries(
        network_time_series[:, network - 1:network, :, :], time,
        ConditionLabels=CONDITION_LABELS,
        SignificantWindows=[anova_windows, *pair_windows],
        SignificanceLabels=["A", *[str(i) for i in range(1, len(condition_pairs) + 1)]],
        SignificanceStyle="line", FigureMode=FIGURE_MODE,
        OutputPath=output_path, FigurePrefix=f"Within_BN{network}",
    )
    statistical_figures.append(figure)
CONDITION_STATS.Figures = statistical_figures
print("A = overall condition ANOVA; numbered lines =", contrast_labels)


# %% 4) CONDITION DIFFERENCES IN SPATIAL ACTIVATION PATTERNS
RUN_SPATIAL_STATISTICS = False  # Set True after the time-series lesson.
if RUN_SPATIAL_STATISTICS:
    SPATIAL_STATS, SPATIAL_FIGURES = BROADNESS_SpatialPatternStatistics(
        BROADNESS, principalcomps=SELECTED_NETWORKS,
        conditionpairs=condition_pairs, conditionnames=CONDITION_LABELS,
        mni_coords=mni_coords, permutations=20 if TUTORIAL_MODE else 5000,
        clusteralpha=0.05, alpha=ALPHA, figuremode="off",
    )
    print("Spatial tests complete. Inspect SPATIAL_STATS.Voxelwise.")
else:
    SPATIAL_STATS = Result(Enabled=False)
    SPATIAL_FIGURES = None
    print("Spatial tests skipped. They use resampling across source maps.")


# %% 5) OPTIONAL ASSOCIATION WITH BEHAVIOUR
# Enter one real behavioural score per participant, in participant_names order.
BEHAVIOUR = None
if BEHAVIOUR is not None:
    behaviour = np.asarray(BEHAVIOUR, dtype=float)
    if behaviour.shape != (n_participants,):
        raise ValueError("BEHAVIOUR must have one score per participant in printed file order")
    correlations, behaviour_p, behaviour_adjusted, behaviour_significant, behaviour_critical = (
        BROADNESS_BehaviourCorrelation(selected_series, behaviour, ALPHA, "Pearson")
    )
    BEHAVIOUR_STATS = Result(Enabled=True, Correlations=correlations, PValues=behaviour_p,
                             AdjustedPValues=behaviour_adjusted,
                             Significant=behaviour_significant, FDRCriticalP=behaviour_critical)
    print("Behavioural association complete. Inspect BEHAVIOUR_STATS.")
else:
    BEHAVIOUR_STATS = Result(Enabled=False)
    print("Behavioural association skipped; no scores were supplied.")


# %% 6) OPTIONAL PHASE-SPACE OR RQA ASSOCIATION WITH BEHAVIOUR
PHASE_BEHAVIOUR_ANALYSIS = "RQA"  # Or "PhaseSpace".
if BEHAVIOUR is not None:
    RQA_BROADNESS, _ = BROADNESS_PhaseSpace_RQA(
        BROADNESS, principalcomps=SELECTED_NETWORKS,
        figuremode="off", video="off",
    )
    if PHASE_BEHAVIOUR_ANALYSIS.lower() == "rqa":
        metric_names = list(RQA_BROADNESS.RQA_metrics[0].columns)
        metric_data = np.stack([table.to_numpy().T for table in RQA_BROADNESS.RQA_metrics], axis=-1)
    elif PHASE_BEHAVIOUR_ANALYSIS.lower() == "phasespace":
        PHASE_STATS = BROADNESS_PhaseSpaceStatistics(
            RQA_BROADNESS, conditionnames=CONDITION_LABELS, statistics="off",
        )
        metric_names = ["Speed", "Dispersion"]
        metric_data = np.stack([getattr(PHASE_STATS.Metrics, name) for name in metric_names], axis=1)
    else:
        raise ValueError("PHASE_BEHAVIOUR_ANALYSIS must be 'RQA' or 'PhaseSpace'")
    r, p, adjusted, sig, critical = BROADNESS_BehaviourCorrelation(metric_data, behaviour, ALPHA)
    PHASE_RQA_BEHAVIOUR = Result(Enabled=True, Analysis=PHASE_BEHAVIOUR_ANALYSIS,
                                 MetricNames=metric_names, Correlations=r, PValues=p,
                                 AdjustedPValues=adjusted, Significant=sig,
                                 FDRCriticalP=critical)
    print("Phase/RQA behavioural association complete.")
else:
    PHASE_RQA_BEHAVIOUR = Result(Enabled=False)
    print("Phase/RQA behavioural association skipped; no scores were supplied.")
