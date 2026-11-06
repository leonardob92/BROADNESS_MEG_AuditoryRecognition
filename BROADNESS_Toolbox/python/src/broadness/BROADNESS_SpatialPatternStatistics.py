"""Paired spatial cluster permutation inference on individual network maps."""

from __future__ import annotations

from types import SimpleNamespace as Result

import numpy as np
from numpy.typing import ArrayLike
from scipy.stats import f as f_dist, t as t_dist

from .BROADNESS_FigureSettings import BROADNESS_FigureSettings
from .BROADNESS_RepeatedMeasuresANOVA import BROADNESS_RepeatedMeasuresANOVA
from ._spatial import (
    SUMMARY_NAMES, cluster_outputs, coordinates, fdr, label_clusters,
    maps_and_components, names_or_default, paired_t, pairs_or_default,
    prepare_maps, repeated_f, summaries, validate_options,
)


def BROADNESS_SpatialPatternStatistics(
    BROADNESS: object, *, principalcomps: ArrayLike = (1, 2),
    conditionpairs: ArrayLike | None = None, conditionnames: list[str] | None = None,
    mni_coords: ArrayLike | None = None, mapnormalization: str = "none",
    omnibus: str = "off", permutations: int = 1000, clusteralpha: float = 0.01,
    alpha: float = 0.05, neighbourdistance: float | None = None,
    randomseed: int = 1, figuremode: str = "off", figurelayout: str = "both",
    OutputPath: str | None = None, figureformats: tuple[str, ...] = ("png",),
    figureprefix: str = "",
) -> tuple[Result, Result]:
    """Test paired condition maps with participant sign flips.

    The maximum cluster mass spans all selected components within each
    condition contrast. A cluster's corrected p-value is repeated over its
    sources; individual sources do not each have that p-value.
    """
    validate_options(permutations, clusteralpha, alpha, randomseed, figuremode)
    if omnibus not in ("on", "off"):
        raise ValueError("omnibus must be 'on' or 'off'")
    figure_settings = BROADNESS_FigureSettings(figuremode, figurelayout, OutputPath,
                                               figureformats, figureprefix,
                                               "SpatialPatternStatistics")
    maps, pcs = maps_and_components(BROADNESS, principalcomps)
    n_sources, _, n_conditions, n_participants = maps.shape
    if n_conditions < 2 or n_participants < 2:
        raise ValueError("at least two conditions and participants are required")
    pairs = pairs_or_default(conditionpairs, n_conditions)
    names = names_or_default(conditionnames, n_conditions, "Condition")
    coords = coordinates(mni_coords, n_sources)
    valid, analysis, neighbours, distance = prepare_maps(maps, pcs, coords, mapnormalization, neighbourdistance)
    raw = maps[valid][:, pcs - 1]
    n_valid, n_pcs, _, _ = analysis.shape
    n_pairs = len(pairs)
    shape = (n_sources, n_pcs, n_pairs)
    mean_difference = np.full(shape, np.nan)
    t_map = np.full(shape, np.nan)
    dz_map = np.full(shape, np.nan)
    p_map = np.full(shape, np.nan)
    cluster_labels = np.zeros(shape, dtype=int)
    cluster_p = np.full(shape, np.nan)
    significant = np.zeros(shape, dtype=bool)
    clusters = np.empty((n_pcs, n_pairs), dtype=object)
    null_mass = np.zeros((permutations, n_pairs))
    threshold = float(t_dist.ppf(1 - clusteralpha / 2, n_participants - 1))
    generator = np.random.default_rng(randomseed)

    for pair_index, (first, second) in enumerate(pairs):
        differences = analysis[:, :, first - 1] - analysis[:, :, second - 1]
        for component in range(n_pcs):
            t, delta, dz = paired_t(differences[:, component, :])
            mean_difference[valid, component, pair_index] = delta
            t_map[valid, component, pair_index] = t
            dz_map[valid, component, pair_index] = dz
            p_map[valid, component, pair_index] = 2 * t_dist.sf(np.abs(t), n_participants - 1)
        signs = generator.choice((-1.0, 1.0), size=(permutations, n_participants))
        for b in range(permutations):
            for component in range(n_pcs):
                t, _, _ = paired_t(differences[:, component, :] * signs[b])
                _, masses, _ = label_clusters(t, threshold, neighbours, True)
                if masses.size:
                    null_mass[b, pair_index] = max(null_mass[b, pair_index], np.max(masses))
        for component in range(n_pcs):
            labels, corrected, sig, rows = cluster_outputs(
                t_map[valid, component, pair_index], threshold, neighbours,
                null_mass[:, pair_index], valid, coords, alpha, True,
            )
            cluster_labels[:, component, pair_index] = labels
            cluster_p[:, component, pair_index] = corrected
            significant[:, component, pair_index] = sig
            clusters[component, pair_index] = rows

    omni = Result(Enabled=omnibus == "on")
    if omni.Enabled:
        df1 = n_conditions - 1
        df2 = (n_participants - 1) * df1
        f_threshold = float(f_dist.ppf(1 - clusteralpha, df1, df2))
        f_map = np.full((n_sources, n_pcs), np.nan)
        f_p = np.full_like(f_map, np.nan)
        f_labels = np.zeros_like(f_map, dtype=int)
        f_cluster_p = np.full_like(f_map, np.nan)
        f_sig = np.zeros_like(f_map, dtype=bool)
        f_clusters = np.empty(n_pcs, dtype=object)
        f_null = np.zeros(permutations)
        for component in range(n_pcs):
            f = repeated_f(analysis[:, component])
            f_map[valid, component] = f
            f_p[valid, component] = f_dist.sf(f, df1, df2)
        # The same condition-label permutation is used for all sources and components.
        for b in range(permutations):
            orders = np.stack([generator.permutation(n_conditions) for _ in range(n_participants)], axis=1)
            for component in range(n_pcs):
                values = analysis[:, component]
                permuted = np.empty_like(values)
                for participant in range(n_participants):
                    permuted[:, :, participant] = values[:, orders[:, participant], participant]
                f = repeated_f(permuted)
                _, masses, _ = label_clusters(f, f_threshold, neighbours, False)
                if masses.size:
                    f_null[b] = max(f_null[b], np.max(masses))
        for component in range(n_pcs):
            labels, corrected, sig, rows = cluster_outputs(
                f_map[valid, component], f_threshold, neighbours, f_null, valid, coords, alpha, False,
            )
            f_labels[:, component], f_cluster_p[:, component], f_sig[:, component] = labels, corrected, sig
            f_clusters[component] = rows
        omni = Result(Enabled=True, FStatistics=f_map, PValues=f_p, ClusterLabels=f_labels,
                      ClusterCorrectedP=f_cluster_p, Significant=f_sig, Clusters=f_clusters,
                      NullMaximumClusterMass=f_null, ClusterFormingThreshold=f_threshold,
                      DegreesFreedom=np.array([df1, df2]),
                      CorrectionFamily="All valid sources and selected components")
    else:
        omni.Reason = "The optional omnibus condition test was switched off."

    summary_values = summaries(raw, coords[valid])
    summary_shape = (6, n_pcs, n_pairs)
    summary_t = np.full(summary_shape, np.nan)
    summary_dz = np.full(summary_shape, np.nan)
    summary_p = np.full(summary_shape, np.nan)
    summary_delta = np.full(summary_shape, np.nan)
    for pair_index, (first, second) in enumerate(pairs):
        differences = summary_values[:, :, first - 1] - summary_values[:, :, second - 1]
        for measure in range(6):
            for component in range(n_pcs):
                values = differences[measure, component]
                values = values[np.isfinite(values)]
                if len(values) < 2:
                    continue
                mean = values.mean()
                sd = values.std(ddof=1)
                t = mean / (sd / np.sqrt(len(values))) if sd > 0 else (np.sign(mean) * np.inf if mean else 0.0)
                summary_t[measure, component, pair_index] = t
                summary_p[measure, component, pair_index] = 2 * t_dist.sf(abs(t), len(values) - 1)
                summary_delta[measure, component, pair_index] = mean
                summary_dz[measure, component, pair_index] = mean / sd if sd > 0 else np.nan
    summary_sig, summary_adj, summary_crit = fdr(summary_p, alpha)
    summary_anova = (BROADNESS_RepeatedMeasuresANOVA(summary_values, alpha)
                     if omnibus == "on" else Result(Enabled=False, Reason="The optional omnibus condition test was switched off."))
    summary = Result(Names=list(SUMMARY_NAMES), Values=summary_values)
    for i, name in enumerate(SUMMARY_NAMES):
        setattr(summary, name, summary_values[i])
    figures = Result(Handles=[], Files=[])
    pairwise = Result(
        ConditionPairs=pairs, ContrastLabels=[f"{names[a - 1]} vs {names[b - 1]}" for a, b in pairs],
        MeanDifference=mean_difference, TStatistics=t_map, CohenDz=dz_map,
        PValues=p_map, ClusterLabels=cluster_labels, ClusterCorrectedP=cluster_p,
        Significant=significant, Clusters=clusters, NullMaximumClusterMass=null_mass,
        ClusterFormingThreshold=threshold, DegreesFreedom=n_participants - 1,
        CorrectionFamily="All valid sources and selected components, separately per condition pair",
    )
    summary_pairwise = Result(
        ConditionPairs=pairs, ContrastLabels=pairwise.ContrastLabels,
        MeanDifference=summary_delta, TStatistics=summary_t, CohenDz=summary_dz,
        PValues=summary_p, AdjustedPValues=summary_adj, Significant=summary_sig,
        FDRCriticalP=summary_crit, FDRMethod="Benjamini-Hochberg",
        FDRFamily="All summary measures, selected components, and condition pairs",
    )
    settings = Result(SelectedComponents=pcs, ConditionPairs=pairs, ConditionNames=names,
                      MapNormalization=mapnormalization, Permutations=permutations,
                      ClusterFormingAlpha=clusteralpha, Alpha=alpha,
                      NeighbourDistance=distance, RandomSeed=randomseed, Omnibus=omnibus,
                      InferenceData="Unthresholded individual spatial patterns",
                      ClusterCorrection="Maximum cluster mass")
    result = Result(MNI_Coordinates=coords, ValidSources=valid,
                    Voxelwise=Result(Pairwise=pairwise, Omnibus=omni),
                    Summary=summary, SummaryTests=Result(Pairwise=summary_pairwise, ANOVA=summary_anova),
                    Settings=settings, Figures=figures)
    if figure_settings.Mode != "off":
        from ._analysis_figures import paired_spatial_figures

        figures = paired_spatial_figures(result, figure_settings)
        result.Figures = figures
    return result, figures


__all__ = ["BROADNESS_SpatialPatternStatistics"]
