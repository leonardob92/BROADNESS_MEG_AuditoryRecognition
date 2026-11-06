"""Cluster source activation profiles across selected brain networks."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path
from types import SimpleNamespace as Result

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from .BROADNESS_FigureSettings import BROADNESS_FigureSettings
from ._spatial import coordinates, field


def _ordered_solution(
    labels: np.ndarray, centroids: np.ndarray, activations: np.ndarray,
    pcs: np.ndarray, support_threshold: float,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    k = len(centroids)
    scale = np.maximum(np.max(np.abs(activations), axis=0), np.finfo(float).eps)
    support = np.empty((k, len(pcs)))
    strength = np.empty_like(support)
    signed_mean = np.empty_like(support)
    counts = np.empty(k, dtype=int)
    first_source = np.empty(k, dtype=int)
    for cluster in range(k):
        mask = labels == cluster
        values = activations[mask]
        support[cluster] = np.mean(values != 0, axis=0)
        strength[cluster] = np.mean(np.abs(values), axis=0) / scale
        signed_mean[cluster] = np.mean(values, axis=0)
        counts[cluster] = mask.sum()
        first_source[cluster] = np.flatnonzero(mask)[0]
    represented = support >= support_threshold
    breadth = represented.sum(axis=1)
    total_support = support.sum(axis=1)
    total_strength = strength.sum(axis=1)
    sort_columns = np.column_stack((breadth, total_support, total_strength, represented, strength, -first_source))
    order = sorted(range(k), key=lambda i: tuple(-v for v in sort_columns[i]))
    mapping = np.empty(k, dtype=int)
    mapping[order] = np.arange(1, k + 1)
    relabeled = mapping[labels]
    combinations = [" + ".join(f"BN{pc}" for pc, flag in zip(pcs, represented[i]) if flag) or "None/background" for i in order]
    summary = pd.DataFrame({
        "Cluster": np.arange(1, k + 1), "OriginalCluster": np.asarray(order) + 1,
        "VoxelCount": counts[order], "RepresentedNetworks": breadth[order],
        "TotalSupport": total_support[order], "TotalNormalizedStrength": total_strength[order],
        "NetworkCombination": combinations,
    })
    for column, pc in enumerate(pcs):
        summary[f"BN{pc}_Represented"] = represented[order, column]
        summary[f"BN{pc}_Support"] = support[order, column]
        summary[f"BN{pc}_NormalizedStrength"] = strength[order, column]
        summary[f"BN{pc}_SignedMean"] = signed_mean[order, column]
    return relabeled, centroids[order], summary


def BROADNESS_SpatialActivationClustering(
    BROADNESS: object, *, principalcomps: ArrayLike = (1, 2),
    nclusters: ArrayLike = tuple(range(2, 21)), evalclusters: int = 10,
    thresh: float | None = None, representation_threshold: float = 0.5,
    scatterplots: str | None = None, OutputPath: str | Path | None = None,
    outpath: str | Path | None = None, mni_coords: ArrayLike | None = None,
    brainmarkersize: float = 8, braincolorintensity: float = 0.75,
    figuremode: str = "off", figurelayout: str = "individual",
    figureoutpath: str | Path | None = None,
    figureformats: tuple[str, ...] = ("png",), figureprefix: str = "",
    kmeans_restarts: int = 10, silhouette_sample_size: int = 500,
) -> tuple[Result, Result]:
    """Cluster thresholded source loadings, ordering clusters by network breadth.

    Python uses k-means++ with ten restarts by default and samples up to 500
    sources for silhouette selection. These practical defaults reduce the
    large cost of MATLAB's 100 restarts and full silhouette calculation;
    increase both options when reproducing an analysis closely.
    """
    if scatterplots not in (None, "all"):
        raise ValueError("scatterplots must be None or 'all'")
    if not brainmarkersize > 0 or not 0 < braincolorintensity <= 1:
        raise ValueError("brainmarkersize and braincolorintensity must be positive")
    figure_settings = BROADNESS_FigureSettings(
        figuremode, figurelayout,
        OutputPath if OutputPath is not None else (outpath if outpath is not None else figureoutpath),
        figureformats, figureprefix, "SpatialActivationClustering",
    )
    raw = np.asarray(field(BROADNESS, "ActivationPatterns_BrainNetworks"), dtype=float)
    original = np.asarray(field(BROADNESS, "OriginalData"))
    if raw.ndim != 2 or original.shape[0] != raw.shape[0] or not np.all(np.isfinite(raw)):
        raise ValueError("activation patterns must be finite sources × components and match OriginalData")
    n_sources, n_components = raw.shape
    pcs = np.asarray(principalcomps)
    if pcs.ndim != 1 or not pcs.size or not np.all(np.isfinite(pcs)) or np.any(pcs != np.floor(pcs)) or np.any((pcs < 1) | (pcs > n_components)) or np.unique(pcs).size != pcs.size:
        raise ValueError("principalcomps must contain unique available one-based components")
    pcs = pcs.astype(int)
    ks = np.asarray(nclusters)
    if ks.ndim != 1 or not ks.size or not np.all(np.isfinite(ks)) or np.any(ks != np.floor(ks)) or np.any((ks < 2) | (ks >= n_sources)) or np.unique(ks).size != ks.size:
        raise ValueError("nclusters must contain unique integers between 2 and n_sources-1")
    ks = ks.astype(int)
    if not 0 < representation_threshold <= 1:
        raise ValueError("representation_threshold must be in (0,1]")
    if not isinstance(evalclusters, int) or evalclusters < 1 or kmeans_restarts < 1 or silhouette_sample_size < 2:
        raise ValueError("evalclusters, kmeans_restarts, and silhouette_sample_size must be positive")
    coords = coordinates(mni_coords, n_sources)
    output = OutputPath if OutputPath is not None else outpath
    if output is not None:
        try:
            bundled_coords = coordinates(None, n_sources)
        except ValueError:
            raise ValueError("NIfTI export requires sources in the bundled MNI voxel order") from None
        if not np.allclose(coords, bundled_coords):
            raise ValueError("NIfTI export requires sources in the bundled MNI voxel order")
    selected = raw[:, pcs - 1]
    if thresh is None:
        thresholds = np.mean(np.abs(selected), axis=0) + np.std(np.abs(selected), axis=0, ddof=1)
    else:
        thresholds = np.broadcast_to(np.asarray(thresh, dtype=float), (pcs.size,))
        if np.any(~np.isfinite(thresholds)) or np.any(thresholds < 0):
            raise ValueError("thresh must be a nonnegative scalar or one value per component")
    activations = np.where(np.abs(selected) > thresholds, selected, 0)
    sd = activations.std(axis=0, ddof=1)
    if np.any(sd == 0):
        raise ValueError("a thresholded component has zero variance; reduce thresh or select another component")
    standardized = (activations - activations.mean(axis=0)) / sd
    labels_by_k: dict[int, np.ndarray] = {}
    centroids_by_k: dict[int, np.ndarray] = {}
    summaries_by_k: dict[int, pd.DataFrame] = {}
    inertia_by_k: dict[int, float] = {}
    for k in ks:
        solution = KMeans(n_clusters=int(k), n_init=kmeans_restarts, random_state=42).fit(standardized)
        labels, centroids, summary = _ordered_solution(solution.labels_, solution.cluster_centers_, activations, pcs, representation_threshold)
        labels_by_k[int(k)] = labels
        centroids_by_k[int(k)] = centroids
        summaries_by_k[int(k)] = summary
        inertia_by_k[int(k)] = solution.inertia_
    if len(ks) == 1:
        optimal_k = int(ks[0])
    else:
        winners = []
        for repetition in range(1, evalclusters + 1):
            scores = []
            for k in ks:
                labels = KMeans(n_clusters=int(k), n_init=kmeans_restarts, random_state=repetition).fit_predict(standardized)
                sample = min(n_sources, silhouette_sample_size)
                score = silhouette_score(standardized, labels, sample_size=sample if sample < n_sources else None, random_state=repetition)
                scores.append(score)
            winners.append(int(ks[np.argmax(scores)]))
        counts = np.bincount(winners, minlength=int(max(ks)) + 1)
        optimal_k = int(np.flatnonzero(counts == counts.max())[0])
    labels = labels_by_k[optimal_k]
    pc_names = [f"PC{pc}" for pc in pcs]
    cluster_info = []
    points = []
    minmax = []
    for cluster in range(1, optimal_k + 1):
        chosen = labels == cluster
        values = activations[chosen]
        cluster_info.append(pd.DataFrame(np.column_stack((np.arange(1, chosen.sum() + 1), coords[chosen], values)), columns=["Voxel", "X", "Y", "Z", *pc_names]))
        points.append(pd.DataFrame(values, columns=pc_names))
        row = {"Cluster": cluster}
        for column, name in enumerate(pc_names):
            row[f"{name}_min"] = values[:, column].min()
            row[f"{name}_max"] = values[:, column].max()
        minmax.append(row)
    names = [f"Nclusters_{k}" for k in ks]
    result = Result(
        idx=pd.DataFrame({name: labels_by_k[int(k)] for name, k in zip(names, ks)}),
        SUM=pd.DataFrame([{name: inertia_by_k[int(k)] for name, k in zip(names, ks)}]),
        Centroids=pd.DataFrame([{name: centroids_by_k[int(k)] for name, k in zip(names, ks)}]),
        ClusterSummaryAll=pd.DataFrame([{name: summaries_by_k[int(k)] for name, k in zip(names, ks)}]),
        RepresentationThreshold=representation_threshold,
        ActivationThresholds=pd.DataFrame([{f"BN{pc}": value for pc, value in zip(pcs, thresholds)}]),
        SelectedComponents=pcs, optimalK=optimal_k,
        ClusterSummary=summaries_by_k[optimal_k], Clusters_info=cluster_info,
        ClusterPoints_PC=points, ClusterMinMax_PC=pd.DataFrame(minmax),
        Figures=Result(Handles=[], Files=[]),
    )
    if output is not None:
        import nibabel as nib

        template = nib.load(str(files("broadness").joinpath("assets/MNI152_8mm_brain_diy.nii.gz")))
        ids = np.asarray(template.dataobj).astype(int)
        valid_ids = (ids >= 1) & (ids <= n_sources)
        destination = Path(output) / "BROADNESS_nifti"
        destination.mkdir(parents=True, exist_ok=True)
        result.NIfTIFiles = []
        for cluster in range(1, optimal_k + 1):
            volume = np.zeros(ids.shape, dtype=np.uint8)
            volume[valid_ids] = (labels[ids[valid_ids] - 1] == cluster).astype(np.uint8)
            path = destination / f"SpatialActivationClustering_OptimalK_{optimal_k}_Cluster_{cluster}.nii.gz"
            nib.save(nib.Nifti1Image(volume, template.affine, template.header), str(path))
            result.NIfTIFiles.append(path)
    if figure_settings.Mode != "off":
        from ._analysis_figures import clustering_figures

        result.Figures = clustering_figures(result, coords, figure_settings, scatterplots,
                                            brainmarkersize, braincolorintensity)
    return result, result.Figures


__all__ = ["BROADNESS_SpatialActivationClustering"]
