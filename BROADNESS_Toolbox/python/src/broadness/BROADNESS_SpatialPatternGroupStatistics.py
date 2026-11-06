"""Independent-group spatial maps with Welch tests and wild bootstrap."""

from __future__ import annotations

from types import SimpleNamespace as Result

import numpy as np
from numpy.typing import ArrayLike
from scipy.stats import f as f_dist, t as t_dist

from .BROADNESS_FigureSettings import BROADNESS_FigureSettings
from ._spatial import (
    SUMMARY_NAMES, cluster_outputs, coordinates, fdr, label_clusters,
    maps_and_components, names_or_default, pairs_or_default, prepare_maps,
    summaries, validate_options, welch_f, welch_t,
)


def BROADNESS_SpatialPatternGroupStatistics(
    BROADNESS: object, group_labels: ArrayLike, *,
    principalcomps: ArrayLike = (1, 2), conditionnames: list[str] | None = None,
    groupnames: list[str] | None = None, mni_coords: ArrayLike | None = None,
    mapnormalization: str = "none", permutations: int = 1000,
    clusteralpha: float = 0.01, alpha: float = 0.05,
    neighbourdistance: float | None = None, randomseed: int = 1,
    figuremode: str = "off", figurelayout: str = "both",
    OutputPath: str | None = None, figureformats: tuple[str, ...] = ("png",),
    figureprefix: str = "",
) -> tuple[Result, Result]:
    """Compare independent participants within a common decomposition.

    The group-centered Rademacher wild bootstrap is approximate. Each
    condition and contrast has a separate maximum-mass family spanning all
    selected sources and components. Use many more than the small resample
    counts shown in tests for final scientific inference.
    """
    validate_options(permutations, clusteralpha, alpha, randomseed, figuremode)
    figure_settings = BROADNESS_FigureSettings(figuremode, figurelayout, OutputPath,
                                               figureformats, figureprefix,
                                               "SpatialPatternGroupStatistics")
    maps, pcs = maps_and_components(BROADNESS, principalcomps)
    n_sources, _, n_conditions, n_participants = maps.shape
    labels = np.asarray(group_labels)
    if labels.ndim != 1 or labels.size != n_participants or not np.all(np.isfinite(labels)) or np.any(labels != np.floor(labels)) or np.any(labels < 1):
        raise ValueError("group_labels must contain one positive integer per participant")
    labels = labels.astype(int)
    n_groups = int(labels.max())
    if n_groups < 2 or not np.array_equal(np.unique(labels), np.arange(1, n_groups + 1)):
        raise ValueError("group labels must be consecutive from 1")
    members = [np.flatnonzero(labels == i) for i in range(1, n_groups + 1)]
    if any(len(indices) < 3 for indices in members):
        raise ValueError("each group needs at least three participants")
    condition_names = names_or_default(conditionnames, n_conditions, "Condition")
    group_names = names_or_default(groupnames, n_groups, "Group")
    pairs = pairs_or_default(None, n_groups)
    coords = coordinates(mni_coords, n_sources)
    valid, analysis, neighbours, distance = prepare_maps(maps, pcs, coords, mapnormalization, neighbourdistance)
    raw = maps[valid][:, pcs - 1]
    n_pcs, n_pairs = len(pcs), len(pairs)
    shape = (n_sources, n_pcs, n_conditions, n_pairs)
    mean_difference = np.full(shape, np.nan)
    t_map = np.full(shape, np.nan)
    g_map = np.full(shape, np.nan)
    p_map = np.full(shape, np.nan)
    cluster_labels = np.zeros(shape, dtype=int)
    cluster_p = np.full(shape, np.nan)
    significant = np.zeros(shape, dtype=bool)
    clusters = np.empty((n_pcs, n_conditions, n_pairs), dtype=object)
    null_mass = np.zeros((permutations, n_conditions, n_pairs))
    omni_enabled = n_groups > 2
    omni_shape = (n_sources, n_pcs, n_conditions)
    f_map = np.full(omni_shape, np.nan)
    f_p = np.full(omni_shape, np.nan)
    f_labels = np.zeros(omni_shape, dtype=int)
    f_cluster_p = np.full(omni_shape, np.nan)
    f_sig = np.zeros(omni_shape, dtype=bool)
    f_clusters = np.empty((n_pcs, n_conditions), dtype=object)
    f_clusters.fill(None)
    f_null = np.zeros((permutations, n_conditions))
    generator = np.random.default_rng(randomseed)

    for condition in range(n_conditions):
        for pair_index, (first, second) in enumerate(pairs):
            x = analysis[:, :, condition, :][:, :, members[first - 1]]
            y = analysis[:, :, condition, :][:, :, members[second - 1]]
            t, df, delta, g = welch_t(x, y)
            t_map[valid, :, condition, pair_index] = t
            mean_difference[valid, :, condition, pair_index] = delta
            g_map[valid, :, condition, pair_index] = g
            p_map[valid, :, condition, pair_index] = 2 * t_dist.sf(np.abs(t), df)
            null_x = (x - x.mean(axis=-1, keepdims=True)) * np.sqrt(x.shape[-1] / (x.shape[-1] - 1))
            null_y = (y - y.mean(axis=-1, keepdims=True)) * np.sqrt(y.shape[-1] / (y.shape[-1] - 1))
            for b in range(permutations):
                e1 = generator.choice((-1.0, 1.0), size=x.shape[-1])
                e2 = generator.choice((-1.0, 1.0), size=y.shape[-1])
                bt, bdf, _, _ = welch_t(null_x * e1, null_y * e2)
                for component in range(n_pcs):
                    threshold = t_dist.ppf(1 - clusteralpha / 2, bdf[:, component])
                    _, masses, _ = label_clusters(bt[:, component], threshold, neighbours, True)
                    if masses.size:
                        null_mass[b, condition, pair_index] = max(null_mass[b, condition, pair_index], np.max(masses))
            for component in range(n_pcs):
                threshold = t_dist.ppf(1 - clusteralpha / 2, df[:, component])
                lab, corrected, sig, rows = cluster_outputs(
                    t[:, component], threshold, neighbours, null_mass[:, condition, pair_index], valid, coords, alpha, True,
                )
                cluster_labels[:, component, condition, pair_index] = lab
                cluster_p[:, component, condition, pair_index] = corrected
                significant[:, component, condition, pair_index] = sig
                clusters[component, condition, pair_index] = rows
        if omni_enabled:
            groups = [analysis[:, :, condition, :][:, :, ids] for ids in members]
            observed_f, observed_df = welch_f(groups)
            f_map[valid, :, condition] = observed_f
            f_p[valid, :, condition] = f_dist.sf(observed_f, n_groups - 1, observed_df)
            residuals = [(group - group.mean(axis=-1, keepdims=True)) * np.sqrt(group.shape[-1] / (group.shape[-1] - 1)) for group in groups]
            for b in range(permutations):
                multipliers = [generator.choice((-1.0, 1.0), size=group.shape[-1]) for group in groups]
                boot_f, boot_df = welch_f([group * e for group, e in zip(residuals, multipliers)])
                for component in range(n_pcs):
                    threshold = f_dist.ppf(1 - clusteralpha, n_groups - 1, boot_df[:, component])
                    _, masses, _ = label_clusters(boot_f[:, component], threshold, neighbours, False)
                    if masses.size:
                        f_null[b, condition] = max(f_null[b, condition], np.max(masses))
            for component in range(n_pcs):
                threshold = f_dist.ppf(1 - clusteralpha, n_groups - 1, observed_df[:, component])
                lab, corrected, sig, rows = cluster_outputs(
                    observed_f[:, component], threshold, neighbours, f_null[:, condition], valid, coords, alpha, False,
                )
                f_labels[:, component, condition] = lab
                f_cluster_p[:, component, condition] = corrected
                f_sig[:, component, condition] = sig
                f_clusters[component, condition] = rows

    summary_values = summaries(raw, coords[valid])
    summary_shape = (6, n_pcs, n_conditions, n_pairs)
    summary_t = np.full(summary_shape, np.nan)
    summary_g = np.full(summary_shape, np.nan)
    summary_p = np.full(summary_shape, np.nan)
    summary_delta = np.full(summary_shape, np.nan)
    summary_f = np.full((6, n_pcs, n_conditions), np.nan)
    summary_f_p = np.full_like(summary_f, np.nan)
    for measure in range(6):
        for component in range(n_pcs):
            for condition in range(n_conditions):
                values = summary_values[measure, component, condition]
                groups = [values[ids][np.isfinite(values[ids])] for ids in members]
                if omni_enabled and all(len(group) >= 2 and group.var(ddof=1) > 0 for group in groups):
                    f, df = welch_f([group[None] for group in groups])
                    summary_f[measure, component, condition] = f[0]
                    summary_f_p[measure, component, condition] = f_dist.sf(f[0], n_groups - 1, df[0])
                for pair_index, (first, second) in enumerate(pairs):
                    x, y = groups[first - 1], groups[second - 1]
                    if len(x) < 2 or len(y) < 2 or x.var(ddof=1) == 0 or y.var(ddof=1) == 0:
                        continue
                    t, df, delta, g = welch_t(x[None], y[None])
                    summary_t[measure, component, condition, pair_index] = t[0]
                    summary_g[measure, component, condition, pair_index] = g[0]
                    summary_p[measure, component, condition, pair_index] = 2 * t_dist.sf(abs(t[0]), df[0])
                    summary_delta[measure, component, condition, pair_index] = delta[0]
    summary_sig, summary_adj, summary_crit = fdr(summary_p, alpha)
    overall_sig, overall_adj, overall_crit = fdr(summary_f_p, alpha)
    summary = Result(Names=list(SUMMARY_NAMES), Values=summary_values)
    for i, name in enumerate(SUMMARY_NAMES):
        setattr(summary, name, summary_values[i])
    pairwise = Result(
        GroupPairs=pairs, MeanDifference=mean_difference, TStatistics=t_map,
        HedgesG=g_map, PValues=p_map, ClusterLabels=cluster_labels,
        ClusterCorrectedP=cluster_p, Significant=significant, Clusters=clusters,
        NullMaximumClusterMass=null_mass,
        CorrectionFamily="Sources and selected components, separately for each condition and group pair",
    )
    omni = Result(
        Enabled=omni_enabled, FStatistics=f_map, PValues=f_p, ClusterLabels=f_labels,
        ClusterCorrectedP=f_cluster_p, Significant=f_sig, Clusters=f_clusters,
        NullMaximumClusterMass=f_null,
        CorrectionFamily="Sources and selected components, separately for each condition",
    )
    summary_pairwise = Result(
        GroupPairs=pairs, MeanDifference=summary_delta, TStatistics=summary_t,
        HedgesG=summary_g, PValues=summary_p, AdjustedPValues=summary_adj,
        Significant=summary_sig, FDRCriticalP=summary_crit, FDRMethod="Benjamini-Hochberg",
        FDRFamily="Measures, selected components, conditions and group pairs",
    )
    summary_overall = Result(
        Enabled=omni_enabled, FStatistics=summary_f, PValues=summary_f_p,
        AdjustedPValues=overall_adj, Significant=overall_sig,
        FDRCriticalP=overall_crit, FDRMethod="Benjamini-Hochberg",
        FDRFamily="Measures, selected components and conditions",
    )
    figures = Result(Handles=[], Files=[])
    settings = Result(
        SelectedComponents=pcs, GroupLabels=labels, GroupNames=group_names, GroupPairs=pairs,
        ConditionNames=condition_names, MapNormalization=mapnormalization,
        Permutations=permutations, ClusterFormingAlpha=clusteralpha, Alpha=alpha,
        NeighbourDistance=distance, RandomSeed=randomseed,
        ResamplingMethod="Group-centered Rademacher wild bootstrap, studentized Welch statistics",
        InferenceData="Unthresholded individual spatial maps",
    )
    result = Result(MNI_Coordinates=coords, ValidSources=valid,
                    Voxelwise=Result(Pairwise=pairwise, Omnibus=omni),
                    Summary=summary,
                    SummaryTests=Result(Pairwise=summary_pairwise, OverallANOVA=summary_overall),
                    Settings=settings, Figures=figures)
    if figure_settings.Mode != "off":
        from ._analysis_figures import group_spatial_figures

        figures = group_spatial_figures(result, figure_settings)
        result.Figures = figures
    return result, figures


__all__ = ["BROADNESS_SpatialPatternGroupStatistics"]
