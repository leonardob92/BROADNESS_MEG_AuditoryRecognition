"""Figure construction for phase-space and spatial analysis results."""

from __future__ import annotations

from math import ceil, sqrt
from types import SimpleNamespace as Result

import numpy as np
from scipy.spatial.distance import pdist, squareform

from .BROADNESS_FinalizeFigure import BROADNESS_FinalizeFigure


def _finish(fig: object, settings: object, figures: Result, name: str) -> None:
    figures.Files.extend(BROADNESS_FinalizeFigure(fig, settings, name))
    if settings.Show:
        figures.Handles.append(fig)


def _phase_axis(ax: object, coords: np.ndarray, time: np.ndarray, condition: int) -> None:
    if coords.shape[1] == 2:
        points = ax.scatter(coords[:, 0], coords[:, 1], c=time, cmap="jet", s=24)
        ax.set_xlabel("Brain network 1")
        ax.set_ylabel("Brain network 2")
    elif coords.shape[1] == 3:
        points = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=time, cmap="jet", s=24)
        ax.set_xlabel("Brain network 1")
        ax.set_ylabel("Brain network 2")
        ax.set_zlabel("Brain network 3")
    else:
        ax.text(0.5, 0.5, "Phase-space display requires 2 or 3 components",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return
    ax.figure.colorbar(points, ax=ax, fraction=0.04, pad=0.08, label="Time (s)")
    ax.set_title(f"Phase space — Condition {condition}")
    ax.grid(True, alpha=0.2)


def _matrix_axis(ax: object, matrix: np.ndarray, time: np.ndarray, title: str, binary: bool) -> None:
    picture = ax.imshow(matrix, origin="lower", aspect="equal",
                        extent=(time[0], time[-1], time[0], time[-1]),
                        cmap="gray" if binary else "viridis_r",
                        vmin=0 if binary else None, vmax=1 if binary else None)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Time (s)")
    ax.set_title(title)
    ax.figure.colorbar(picture, ax=ax, fraction=0.04, pad=0.04)


def rqa_figures(result: Result, threshold: float, theiler_window: int | None, settings: object) -> Result:
    from matplotlib import pyplot as plt

    figures = Result(Handles=[], Files=[])
    time = result.PhaseSpace.Time
    n_conditions = result.PhaseSpace.nConditions
    indices = np.arange(len(time))
    eligible = np.ones((len(time), len(time)), dtype=bool)
    if theiler_window is not None:
        eligible = np.abs(indices[:, None] - indices[None, :]) > theiler_window
    for condition in range(n_conditions):
        coords = result.PhaseSpace.MeanCoordinates[condition, 0]
        distance = squareform(pdist(coords))
        recurrence = ((distance < distance.max() * threshold) & eligible).astype(int)
        projection = "3d" if coords.shape[1] == 3 else None
        if settings.MakeIndividual:
            for kind, name in (("phase", "PhaseSpace"), ("distance", "DistanceMatrix"),
                               ("recurrence", "RecurrenceThresholded")):
                fig = plt.figure(figsize=(6.2, 5.3), constrained_layout=True)
                ax = fig.add_subplot(111, projection=projection if kind == "phase" else None)
                if kind == "phase":
                    _phase_axis(ax, coords, time, condition + 1)
                elif kind == "distance":
                    _matrix_axis(ax, distance, time, f"Condition {condition + 1} — Distance matrix", False)
                else:
                    _matrix_axis(ax, recurrence, time, f"Condition {condition + 1} — Recurrence plot", True)
                _finish(fig, settings, figures, f"Condition_{condition + 1:02d}_{name}")
        if settings.MakeSummary:
            fig = plt.figure(figsize=(15, 4.8), constrained_layout=True)
            fig.suptitle(f"Phase-Space RQA — Condition {condition + 1}")
            ax1 = fig.add_subplot(131, projection=projection)
            ax2 = fig.add_subplot(132)
            ax3 = fig.add_subplot(133)
            _phase_axis(ax1, coords, time, condition + 1)
            _matrix_axis(ax2, distance, time, "Distance matrix", False)
            _matrix_axis(ax3, recurrence, time, "Recurrence plot", True)
            _finish(fig, settings, figures, f"Condition_{condition + 1:02d}_Summary")
    return figures


def rqa_video(result: Result) -> tuple[list[object], list[object]]:
    from matplotlib import pyplot as plt
    from matplotlib.animation import FuncAnimation

    animations = []
    handles = []
    time = result.PhaseSpace.Time
    for condition in range(result.PhaseSpace.nConditions):
        coords = result.PhaseSpace.MeanCoordinates[condition, 0]
        if coords.shape[1] not in (2, 3):
            continue
        fig = plt.figure(figsize=(6.2, 5.3))
        ax = fig.add_subplot(111, projection="3d" if coords.shape[1] == 3 else None)
        ax.set_xlim(coords[:, 0].min(), coords[:, 0].max())
        ax.set_ylim(coords[:, 1].min(), coords[:, 1].max())
        if coords.shape[1] == 3:
            ax.set_zlim(coords[:, 2].min(), coords[:, 2].max())
        ax.set_xlabel("Brain network 1")
        ax.set_ylabel("Brain network 2")
        ax.set_title(f"Phase space — Condition {condition + 1}")
        dots = ax.scatter([], [], [], s=24) if coords.shape[1] == 3 else ax.scatter([], [], s=24)

        def update(frame, points=coords, artist=dots, axis=ax, condition_index=condition):
            if points.shape[1] == 3:
                artist._offsets3d = (points[:frame + 1, 0], points[:frame + 1, 1], points[:frame + 1, 2])
            else:
                artist.set_offsets(points[:frame + 1])
            axis.set_title(f"Phase space — Condition {condition_index + 1}, t={time[frame]:.3f} s")
            return (artist,)

        animation = FuncAnimation(fig, update, frames=len(time), interval=40, blit=False)
        fig.show()
        animations.append(animation)
        handles.append(fig)
    return animations, handles


def _brain_surface(axis: object) -> None:
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    from .BROADNESS_Visualizer import _mni_brain_surface

    vertices, faces = _mni_brain_surface()
    surface = Poly3DCollection(vertices[faces], facecolor=(0.56, 0.61, 0.66, 0.20),
                               edgecolor="none", linewidth=0)
    surface.set_rasterized(True)
    axis.add_collection3d(surface)
    center = (vertices.min(axis=0) + vertices.max(axis=0)) / 2
    radius = max((vertices.max(axis=0) - vertices.min(axis=0)) / 2)
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1, 1, 1), zoom=1.55)


def _map_axes(fig: object, coords: np.ndarray, valid: np.ndarray, stat: np.ndarray,
              significant: np.ndarray, title: str, two_sided: bool) -> None:
    first = fig.add_subplot(121, projection="3d")
    second = fig.add_subplot(122, projection="3d")
    _brain_surface(first)
    _brain_surface(second)
    xyz = coords[valid]
    values = stat[valid]
    limit = np.nanmax(np.abs(values)) if np.any(np.isfinite(values)) else 1
    if not np.isfinite(limit) or limit == 0:
        limit = 1
    dots = first.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=values, s=12,
                         cmap="coolwarm" if two_sided else "viridis",
                         vmin=-limit if two_sided else 0, vmax=limit, alpha=0.8)
    fig.colorbar(dots, ax=first, shrink=0.6, pad=0.04)
    first.set_title("t-statistic" if two_sided else "F-statistic")
    positive = valid & significant & (stat > 0)
    negative = valid & significant & (stat < 0)
    for mask, color in ((positive, (0.76, 0.16, 0.14)), (negative, (0.14, 0.32, 0.72))):
        if np.any(mask) and (two_sided or mask is positive):
            second.scatter(coords[mask, 0], coords[mask, 1], coords[mask, 2],
                           s=22, c=[color])
    if not np.any(significant):
        second.text2D(0.5, 0.05, "No cluster survived correction",
                      transform=second.transAxes, ha="center")
    second.set_title("Cluster-corrected differences")
    for ax in (first, second):
        ax.set_axis_off()
        ax.view_init(elev=10, azim=-90)
    fig.suptitle(title)


def _summary_axes(fig: object, values: np.ndarray, names: list[str], labels: list[str], grouped: bool, members: list[np.ndarray] | None = None) -> None:
    # values: measures × conditions/groups × participants, or measures × participants
    from matplotlib import pyplot as plt

    if values.ndim == 2:
        values = values[:, None, :]
    for measure, name in enumerate(names):
        ax = fig.add_subplot(2, 3, measure + 1)
        for category in range(len(labels)):
            observations = values[measure, category] if members is None else values[measure, 0, members[category]]
            observations = observations[np.isfinite(observations)]
            ax.scatter(np.full(len(observations), category + 1), observations,
                       s=13, alpha=0.45, color=plt.get_cmap("tab10")(category % 10))
            if len(observations):
                se = observations.std(ddof=1) / np.sqrt(len(observations)) if len(observations) > 1 else 0
                ax.errorbar(category + 1, observations.mean(), yerr=se,
                            fmt="ko", mfc="white", capsize=3)
        ax.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=25, ha="right")
        ax.set_xlim(0.5, len(labels) + 0.5)
        ax.set_title(name)
        ax.grid(True, alpha=0.2)


def paired_spatial_figures(result: Result, settings: object) -> Result:
    from matplotlib import pyplot as plt

    figures = Result(Handles=[], Files=[])
    pcs = result.Settings.SelectedComponents
    pairs = result.Settings.ConditionPairs
    valid = result.ValidSources
    coords = result.MNI_Coordinates
    if settings.MakeIndividual:
        for component, pc in enumerate(pcs):
            for pair_index, (first, second) in enumerate(pairs):
                fig = plt.figure(figsize=(11.5, 5), constrained_layout=True)
                _map_axes(fig, coords, valid,
                          result.Voxelwise.Pairwise.TStatistics[:, component, pair_index],
                          result.Voxelwise.Pairwise.Significant[:, component, pair_index],
                          f"Brain Network {pc}: {result.Voxelwise.Pairwise.ContrastLabels[pair_index]}", True)
                _finish(fig, settings, figures, f"SpatialStatistics_BN{pc}_Condition{first}_vs_{second}")
            if result.Voxelwise.Omnibus.Enabled:
                fig = plt.figure(figsize=(11.5, 5), constrained_layout=True)
                _map_axes(fig, coords, valid, result.Voxelwise.Omnibus.FStatistics[:, component],
                          result.Voxelwise.Omnibus.Significant[:, component],
                          f"Brain Network {pc}: overall condition effect", False)
                _finish(fig, settings, figures, f"SpatialStatistics_BN{pc}_Omnibus")
    if settings.MakeSummary:
        for component, pc in enumerate(pcs):
            fig = plt.figure(figsize=(11.5, 7), constrained_layout=True)
            _summary_axes(fig, result.Summary.Values[:, component], result.Summary.Names,
                          result.Settings.ConditionNames, False)
            fig.suptitle(f"Spatial summaries — Brain Network {pc}")
            _finish(fig, settings, figures, f"SpatialSummaries_BN{pc}")
    return figures


def group_spatial_figures(result: Result, settings: object) -> Result:
    from matplotlib import pyplot as plt

    figures = Result(Handles=[], Files=[])
    pcs = result.Settings.SelectedComponents
    coords, valid = result.MNI_Coordinates, result.ValidSources
    pairs = result.Settings.GroupPairs
    members = [np.flatnonzero(result.Settings.GroupLabels == i) for i in range(1, len(result.Settings.GroupNames) + 1)]
    if settings.MakeIndividual:
        for component, pc in enumerate(pcs):
            for condition, condition_name in enumerate(result.Settings.ConditionNames):
                if result.Voxelwise.Omnibus.Enabled:
                    fig = plt.figure(figsize=(11.5, 5), constrained_layout=True)
                    _map_axes(fig, coords, valid, result.Voxelwise.Omnibus.FStatistics[:, component, condition],
                              result.Voxelwise.Omnibus.Significant[:, component, condition],
                              f"Network {pc}: {condition_name} — overall groups", False)
                    _finish(fig, settings, figures, f"SpatialGroup_BN{pc}_C{condition + 1}_Overall")
                for pair_index, (first, second) in enumerate(pairs):
                    fig = plt.figure(figsize=(11.5, 5), constrained_layout=True)
                    _map_axes(fig, coords, valid,
                              result.Voxelwise.Pairwise.TStatistics[:, component, condition, pair_index],
                              result.Voxelwise.Pairwise.Significant[:, component, condition, pair_index],
                              f"Network {pc}: {condition_name} — {result.Settings.GroupNames[first - 1]} vs {result.Settings.GroupNames[second - 1]}", True)
                    _finish(fig, settings, figures, f"SpatialGroup_BN{pc}_C{condition + 1}_G{first}_vs_G{second}")
    if settings.MakeSummary:
        for component, pc in enumerate(pcs):
            for condition, condition_name in enumerate(result.Settings.ConditionNames):
                fig = plt.figure(figsize=(11.5, 7), constrained_layout=True)
                _summary_axes(fig, result.Summary.Values[:, component, condition],
                              result.Summary.Names, result.Settings.GroupNames, True, members)
                fig.suptitle(f"Spatial summaries — Brain Network {pc}: {condition_name}")
                _finish(fig, settings, figures, f"SpatialGroupSummary_BN{pc}_C{condition + 1}")
    return figures


def _cluster_embedding(ax: object, activations: np.ndarray, labels: np.ndarray, pcs: np.ndarray, k: int) -> None:
    if activations.shape[1] not in (2, 3):
        ax.text(0.5, 0.5, "Cluster embedding requires 2 or 3 components",
                transform=ax.transAxes, ha="center")
        ax.set_axis_off()
        return
    colors = np.asarray([plt_color(i, k) for i in labels - 1])
    if activations.shape[1] == 2:
        ax.scatter(activations[:, 0], activations[:, 1], c=colors, s=10, alpha=0.7)
    else:
        ax.scatter(activations[:, 0], activations[:, 1], activations[:, 2], c=colors, s=10, alpha=0.7)
        ax.set_zlabel(f"PC{pcs[2]}")
    ax.set_xlabel(f"PC{pcs[0]}")
    ax.set_ylabel(f"PC{pcs[1]}")
    ax.set_title(f"Cluster embedding (K={k})")


def plt_color(index: int, total: int) -> tuple[float, ...]:
    from matplotlib import pyplot as plt

    return plt.get_cmap("turbo")((index + 0.5) / max(total, 1))


def clustering_figures(result: Result, coords: np.ndarray, settings: object, scatterplots: str | None,
                       brainmarkersize: float, braincolorintensity: float) -> Result:
    from matplotlib import pyplot as plt

    figures = Result(Handles=[], Files=[])
    ks = [int(name.split("_")[-1]) for name in result.idx.columns]
    inertia = result.SUM.iloc[0].to_numpy(dtype=float)
    pcs = result.SelectedComponents
    activations = np.zeros((len(result.idx), len(pcs)))
    for cluster, frame in enumerate(result.ClusterPoints_PC, start=1):
        activations[result.idx[f"Nclusters_{result.optimalK}"].to_numpy() == cluster] = frame.to_numpy()
    def elbow(ax):
        ax.plot(ks, inertia, "o-", color="black")
        ax.axvline(result.optimalK, color="tab:red", linestyle="--", alpha=0.6)
        ax.set_xlabel("Number of clusters")
        ax.set_ylabel("Within-cluster sum of squares")
        ax.set_title("Clustering elbow")
        ax.grid(True, alpha=0.2)
    if settings.MakeIndividual:
        fig, ax = plt.subplots(figsize=(6, 4.5), constrained_layout=True)
        elbow(ax)
        _finish(fig, settings, figures, "ClusteringElbow")
        shown_ks = ks if scatterplots == "all" else [result.optimalK]
        for k in shown_ks:
            fig = plt.figure(figsize=(6, 5), constrained_layout=True)
            ax = fig.add_subplot(111, projection="3d" if len(pcs) == 3 else None)
            _cluster_embedding(ax, activations, result.idx[f"Nclusters_{k}"].to_numpy(), pcs, k)
            _finish(fig, settings, figures, f"ClusterEmbedding_K_{k:02d}")
        for cluster, frame in enumerate(result.ClusterPoints_PC, start=1):
            values = frame.to_numpy()
            mask = result.idx[f"Nclusters_{result.optimalK}"].to_numpy() == cluster
            cluster_color = np.asarray(plt_color(cluster - 1, result.optimalK)[:3]) * braincolorintensity
            title = (f"Cluster {cluster} — {result.ClusterSummary.NetworkCombination.iloc[cluster - 1]}"
                     f" (n = {mask.sum()} sources)")
            fig = plt.figure(figsize=(8.6, 5.8), constrained_layout=True)
            grid = fig.add_gridspec(2, 2, height_ratios=(3.4, 1))
            profile_ax = fig.add_subplot(grid[1, :])
            for column, (elevation, azimuth, label) in enumerate(
                ((18, -60, "Oblique view"), (75, -90, "View from above"))
            ):
                brain_ax = fig.add_subplot(grid[0, column], projection="3d")
                _brain_surface(brain_ax)
                brain_ax.scatter(coords[mask, 0], coords[mask, 1], coords[mask, 2],
                                 s=brainmarkersize, c=[cluster_color], alpha=0.8)
                brain_ax.view_init(elev=elevation, azim=azimuth)
                brain_ax.set_title(label, fontsize=11)
                brain_ax.set_axis_off()
            for column, pc in enumerate(pcs):
                lane = len(pcs) - column
                profile_ax.axhline(lane, color="0.82", linewidth=0.8, zorder=0)
                nonzero = values[:, column][values[:, column] != 0]
                jitter = (np.mod(np.arange(1, len(nonzero) + 1) * 0.61803398875, 1) - 0.5) * 0.32
                profile_ax.scatter(nonzero, lane + jitter, s=12, alpha=0.7, color=cluster_color)
            profile_ax.set_yticks(np.arange(len(pcs), 0, -1), [f"PC{pc}" for pc in pcs])
            profile_ax.set_ylim(0.55, len(pcs) + 0.45)
            profile_ax.set_xlabel("Spatial activation")
            profile_ax.set_title("Activation profile")
            profile_ax.spines[["top", "right", "left"]].set_visible(False)
            profile_ax.tick_params(axis="y", length=0)
            fig.suptitle(title)
            _finish(fig, settings, figures, f"OptimalK_{result.optimalK:02d}_Cluster_{cluster:02d}_ProfileAndBrain")
    if settings.MakeSummary:
        fig = plt.figure(figsize=(11, 4.8), constrained_layout=True)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122, projection="3d" if len(pcs) == 3 else None)
        elbow(ax1)
        _cluster_embedding(ax2, activations, result.idx[f"Nclusters_{result.optimalK}"].to_numpy(), pcs, result.optimalK)
        fig.suptitle(f"Spatial Activation Clustering — Optimal K = {result.optimalK}")
        _finish(fig, settings, figures, f"OptimalK_{result.optimalK:02d}_ClusteringSummary")
    return figures
