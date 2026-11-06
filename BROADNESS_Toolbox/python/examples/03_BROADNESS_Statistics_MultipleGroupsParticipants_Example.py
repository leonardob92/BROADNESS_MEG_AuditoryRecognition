# %% 0) STARTUP AND USER SETTINGS
"""Independent-group BROADNESS statistics on published participant data.

The participant .mat files do not contain research group labels. In TUTORIAL_MODE the
script creates arbitrary demonstration groups so every section can be learned;
these results have no scientific interpretation. Full-data mode requires your
real participant-to-group mapping. See docs/tutorials/03_MULTIPLE_GROUPS.md.
"""

from itertools import combinations
from pathlib import Path
from types import SimpleNamespace as Result

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution, ttest_ind

from broadness import (
    BROADNESS_FDRCorrection,
    BROADNESS_NetworkEstimation,
    BROADNESS_Plot_ActivationTimeseries,
    BROADNESS_SignificantSamplesToWindows,
    BROADNESS_SpatialPatternGroupStatistics,
    BROADNESS_TemporalClusterCorrection_GroupANOVA,
    BROADNESS_TemporalClusterCorrection_Groups,
)
from broadness.example_data import default_data_path, load_published_example

data_path = default_data_path()  # Change to Path("/your/data/folder") if needed.
output_path = default_data_path().parent / "BROADNESS_Toolbox" / "Output" / "Stats"
TUTORIAL_MODE = True
MAX_PARTICIPANTS = 6   # At least 3 in each demonstration group.
ALLOW_LARGE = False    # Explicit opt-in if the selected array exceeds 2 GiB.
PARTICIPANT_FILES = None  # Or a list of exact SUBJ_*.mat names.
GROUP_LABELS_BY_FILE = {}  # For research: {"SUBJ_0002.mat": 1, ...} for every selected file.
GROUP_NAMES = None        # For research: ["Controls", "Patients", ...].
SELECTED_NETWORKS = (1, 2)
CONDITION_LABELS = ["Memorized", "NewT1", "NewT2", "NewT3", "NewT4"]
CORRECTION = "FDR"        # "FDR" or "cluster".
ALPHA = 0.05
CLUSTER_PERMUTATIONS = 100 if TUTORIAL_MODE else 5000
FIGURE_MODE = "show"


# %% 1) LOAD PARTICIPANTS AND ESTIMATE ONE COMMON DECOMPOSITION
example = load_published_example(
    data_path, participants=True, participant_files=PARTICIPANT_FILES,
    max_participants=MAX_PARTICIPANTS, tutorial_mode=TUTORIAL_MODE,
    allow_large=ALLOW_LARGE,
)
data, time, mni_coords = example.data, example.time, example.mni_coords
participant_names = [path.name for path in example.participant_files]
if data.shape[2] != len(CONDITION_LABELS):
    raise ValueError("Update CONDITION_LABELS to match the loaded conditions")
if len(participant_names) < 6:
    raise ValueError("Load at least six participants, allowing three per group")
BROADNESS = BROADNESS_NetworkEstimation(data, time, components_num=max(SELECTED_NETWORKS))
network_time_series = BROADNESS.TimeSeries_BrainNetworks
print("Participant order:", participant_names)
print("Network time-series shape:", network_time_series.shape)


# %% 2) INDEPENDENT GROUP TESTS FOR EACH CONDITION
if GROUP_LABELS_BY_FILE:
    missing_labels = [name for name in participant_names if name not in GROUP_LABELS_BY_FILE]
    if missing_labels:
        raise ValueError(f"Give a group label for every loaded file: {missing_labels}")
    group_labels = np.array([GROUP_LABELS_BY_FILE[name] for name in participant_names], dtype=int)
    if GROUP_NAMES is None:
        raise ValueError("Set GROUP_NAMES for your real group labels")
    group_names = list(GROUP_NAMES)
    demo_groups = False
elif TUTORIAL_MODE:
    midpoint = len(participant_names) // 2
    group_labels = np.array([1] * midpoint + [2] * (len(participant_names) - midpoint))
    group_names = ["Demonstration A", "Demonstration B"]
    demo_groups = True
    print("DEMONSTRATION ONLY: group assignments are arbitrary, not research metadata.")
else:
    raise ValueError("Full-data analysis requires GROUP_LABELS_BY_FILE and GROUP_NAMES")

n_groups = len(group_names)
if n_groups < 2 or not np.array_equal(np.unique(group_labels), np.arange(1, n_groups + 1)):
    raise ValueError("Group labels must be consecutive integers from 1 and match GROUP_NAMES")
group_indices = [np.flatnonzero(group_labels == label) for label in range(1, n_groups + 1)]
if any(len(indices) < 3 for indices in group_indices):
    raise ValueError("Every group needs at least three participants for the spatial section")
group_pairs = tuple(combinations(range(1, n_groups + 1), 2))
selected_series = network_time_series[:, np.array(SELECTED_NETWORKS) - 1, :, :]
n_time, n_networks, n_conditions, _ = selected_series.shape
shape = (n_time, n_networks, n_conditions, len(group_pairs))
p_values, t_statistics, hedges_g = (np.empty(shape) for _ in range(3))

for pair_index, (first, second) in enumerate(group_pairs):
    x = selected_series[:, :, :, group_indices[first - 1]]
    y = selected_series[:, :, :, group_indices[second - 1]]
    test = ttest_ind(x, y, axis=-1, equal_var=False)
    p_values[:, :, :, pair_index] = test.pvalue
    t_statistics[:, :, :, pair_index] = test.statistic
    n1, n2 = x.shape[-1], y.shape[-1]
    pooled_sd = np.sqrt(((n1 - 1) * x.var(axis=-1, ddof=1) +
                         (n2 - 1) * y.var(axis=-1, ddof=1)) / (n1 + n2 - 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        hedges_g[:, :, :, pair_index] = ((x.mean(axis=-1) - y.mean(axis=-1)) /
                                          pooled_sd * (1 - 3 / (4 * (n1 + n2) - 9)))
hedges_g[~np.isfinite(hedges_g)] = np.nan

if CORRECTION.lower() == "fdr":
    significant, adjusted_p, critical_p = BROADNESS_FDRCorrection(p_values, ALPHA)
    cluster_stats = None
elif CORRECTION.lower() == "cluster":
    significant, adjusted_p, critical_p, cluster_stats = BROADNESS_TemporalClusterCorrection_Groups(
        selected_series, group_labels, group_pairs, ALPHA,
        permutations=CLUSTER_PERMUTATIONS, clusteralpha=0.05, time=time,
    )
else:
    raise ValueError("CORRECTION must be 'FDR' or 'cluster'")

if n_groups > 2:
    if CORRECTION.lower() == "cluster":
        overall_significant, overall_p, _, overall_details = BROADNESS_TemporalClusterCorrection_GroupANOVA(
            selected_series, group_labels, ALPHA,
            permutations=CLUSTER_PERMUTATIONS, clusteralpha=0.05, time=time,
        )
    else:
        means = np.stack([selected_series[:, :, :, idx].mean(axis=-1) for idx in group_indices])
        variances = np.stack([selected_series[:, :, :, idx].var(axis=-1, ddof=1) for idx in group_indices])
        counts = np.array([len(idx) for idx in group_indices])[:, None, None, None]
        weights = counts / variances
        total_weight = weights.sum(axis=0)
        weighted_mean = (weights * means).sum(axis=0) / total_weight
        correction_sum = (((1 - weights / total_weight) ** 2) / (counts - 1)).sum(axis=0)
        k = n_groups
        overall_f = ((weights * (means - weighted_mean) ** 2).sum(axis=0) / (k - 1) /
                     (1 + 2 * (k - 2) * correction_sum / (k * k - 1)))
        df2 = (k * k - 1) / (3 * correction_sum)
        overall_raw_p = f_distribution.sf(overall_f, k - 1, df2)
        overall_significant, overall_p, _ = BROADNESS_FDRCorrection(overall_raw_p, ALPHA)
        overall_details = Result(FStatistics=overall_f, PValues=overall_raw_p)
else:
    overall_significant = np.zeros((n_time, n_networks, n_conditions), dtype=bool)
    overall_p = np.full(overall_significant.shape, np.nan)
    overall_details = None

GROUP_STATS = Result(
    Settings=Result(ParticipantNames=participant_names, GroupLabels=group_labels,
                    GroupNames=group_names, ConditionLabels=CONDITION_LABELS,
                    CorrectionMethod=CORRECTION, TutorialMode=TUTORIAL_MODE,
                    DemonstrationGroups=demo_groups),
    BetweenGroupComparisons=Result(GroupPairs=group_pairs, PValues=p_values,
                                   AdjustedPValues=adjusted_p, Significant=significant,
                                   TStatistics=t_statistics, HedgesG=hedges_g,
                                   ClusterResults=cluster_stats,
                                   OverallANOVA=Result(Enabled=n_groups > 2,
                                                       Significant=overall_significant,
                                                       AdjustedPValues=overall_p,
                                                       Details=overall_details)),
)
print("Group tests complete. Inspect GROUP_STATS.BetweenGroupComparisons.")


# %% 3) PLOT SIGNIFICANT TIME WINDOWS
group_members = [list(indices + 1) for indices in group_indices]  # Plot API uses one-based members.
significance_key = [f"{CONDITION_LABELS[c]}: {group_names[a - 1]} vs {group_names[b - 1]}"
                    for c in range(n_conditions) for a, b in group_pairs]
if n_groups > 2:
    significance_key = [f"{name}: overall group ANOVA" for name in CONDITION_LABELS] + significance_key
statistical_figures = []
for network_index, network in enumerate(SELECTED_NETWORKS):
    windows = []
    if n_groups > 2:
        windows.extend(BROADNESS_SignificantSamplesToWindows(
            overall_significant[:, network_index, condition], time
        ) for condition in range(n_conditions))
    windows.extend(BROADNESS_SignificantSamplesToWindows(
        significant[:, network_index, condition, pair_index], time
    ) for condition in range(n_conditions) for pair_index in range(len(group_pairs)))
    figure = BROADNESS_Plot_ActivationTimeseries(
        network_time_series[:, network - 1:network, :, :], time,
        Groups=group_members, GroupLabels=group_names,
        ConditionLabels=CONDITION_LABELS, SignificantWindows=windows,
        SignificanceLabels=[str(i) for i in range(1, len(windows) + 1)],
        SignificanceStyle="line", FigureMode=FIGURE_MODE,
        OutputPath=output_path, FigurePrefix=f"Between_BN{network}",
    )
    statistical_figures.append(figure)
GROUP_STATS.Figures = Result(SignificanceKey=significance_key,
                             GroupComparisons=statistical_figures)
print("Significance-line key:", list(enumerate(significance_key, start=1)))


# %% 4) CONDITION × GROUP INTERACTION IN A PREDEFINED TIME WINDOW
# This mixed model is optional and requires statsmodels (installed by setup).
RUN_INTERACTION = False
SUMMARY_TIME_WINDOW = (0.35, 0.75)  # Choose before inspecting results.
if RUN_INTERACTION:
    if demo_groups:
        raise ValueError(
            "The mixed model needs real GROUP_LABELS_BY_FILE assignments. "
            "Arbitrary demonstration groups are only for sections 2, 3, and 6."
        )
    import statsmodels.formula.api as smf

    window_mask = (time >= SUMMARY_TIME_WINDOW[0]) & (time <= SUMMARY_TIME_WINDOW[1])
    if not window_mask.any():
        raise ValueError("SUMMARY_TIME_WINDOW does not overlap the selected time vector")
    interaction_tables, interaction_models = [], []
    for network in SELECTED_NETWORKS:
        values = network_time_series[window_mask, network - 1].mean(axis=0)
        table = pd.DataFrame([
            {"Response": values[condition, participant], "Participant": str(participant_names[participant]),
             "Condition": str(CONDITION_LABELS[condition]), "Group": str(group_names[group_labels[participant] - 1])}
            for participant in range(len(participant_names)) for condition in range(n_conditions)
        ])
        try:
            model = smf.mixedlm("Response ~ C(Condition) * C(Group)", table,
                                groups=table["Participant"]).fit(reml=False, method="lbfgs")
        except np.linalg.LinAlgError as exc:
            raise RuntimeError(
                "The mixed model is singular for these participants and groups. "
                "Use a suitable larger sample with real group labels, or leave "
                "RUN_INTERACTION = False for the tutorial."
            ) from exc
        interaction_tables.append(table)
        interaction_models.append(model)
        print(f"Network {network} mixed model:\n{model.summary()}")
    GROUP_STATS.ConditionGroupMixedModel = Result(Enabled=True, TimeWindow=SUMMARY_TIME_WINDOW,
                                                  Tables=interaction_tables, Models=interaction_models)
else:
    GROUP_STATS.ConditionGroupMixedModel = Result(Enabled=False)
    print("Mixed model skipped. Set RUN_INTERACTION = True to fit it.")


# %% 5) OPTIONAL BEHAVIOURAL ASSOCIATION
# One real score per participant in participant_names order; None skips.
BEHAVIOUR = None
if BEHAVIOUR is not None:
    import statsmodels.formula.api as smf

    behaviour = np.asarray(BEHAVIOUR, dtype=float)
    if behaviour.shape != (len(participant_names),) or not np.all(np.isfinite(behaviour)):
        raise ValueError("BEHAVIOUR must have one finite score per loaded participant")
    if not GROUP_STATS.ConditionGroupMixedModel.Enabled:
        raise ValueError("Run section 4 with RUN_INTERACTION = True first")
    behaviour_models = []
    for table in GROUP_STATS.ConditionGroupMixedModel.Tables:
        table = table.copy()
        score_by_name = dict(zip(participant_names, behaviour))
        table["Behaviour"] = table["Participant"].map(score_by_name)
        try:
            model = smf.mixedlm("Response ~ C(Condition) * C(Group) + Behaviour",
                                table, groups=table["Participant"]).fit(reml=False, method="lbfgs")
        except np.linalg.LinAlgError as exc:
            raise RuntimeError(
                "The behavioural mixed model is singular. Check sample size, "
                "group assignments, and behaviour variability."
            ) from exc
        behaviour_models.append(model)
    GROUP_STATS.Behaviour = Result(Enabled=True, Models=behaviour_models)
    print("Behavioural mixed models complete. Inspect GROUP_STATS.Behaviour.Models.")
else:
    GROUP_STATS.Behaviour = Result(Enabled=False)
    print("Behavioural models skipped; no scores were supplied.")


# %% 6) GROUP DIFFERENCES IN SPATIAL ACTIVATION PATTERNS
RUN_SPATIAL_STATISTICS = False  # More expensive than the time-series sections.
if RUN_SPATIAL_STATISTICS:
    SPATIAL_GROUP_STATS, SPATIAL_GROUP_FIGURES = BROADNESS_SpatialPatternGroupStatistics(
        BROADNESS, group_labels, principalcomps=SELECTED_NETWORKS,
        groupnames=group_names, conditionnames=CONDITION_LABELS,
        mni_coords=mni_coords, permutations=20 if TUTORIAL_MODE else 5000,
        clusteralpha=0.05, alpha=ALPHA, figuremode="off",
    )
    print("Spatial group tests complete. Inspect SPATIAL_GROUP_STATS.Voxelwise.")
else:
    SPATIAL_GROUP_STATS = Result(Enabled=False)
    SPATIAL_GROUP_FIGURES = None
    print("Spatial group tests skipped.")
