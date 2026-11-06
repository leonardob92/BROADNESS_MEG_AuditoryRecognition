"""Shared numerical helpers for spatial BROADNESS analyses."""

from __future__ import annotations

from importlib.resources import files
from itertools import combinations
from types import SimpleNamespace as Result

import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from scipy.stats import f as f_dist, t as t_dist

from .BROADNESS_FDRCorrection import BROADNESS_FDRCorrection


def field(value: object, name: str) -> object:
    return value[name] if isinstance(value, dict) else getattr(value, name)


def maps_and_components(value: object, principalcomps: object) -> tuple[np.ndarray, np.ndarray]:
    maps = np.asarray(field(value, "ActivationPatterns_BrainNetworks_Individual"), dtype=float)
    if maps.ndim != 4 or not maps.size:
        raise ValueError("individual maps must be sources × components × conditions × participants")
    pcs = np.asarray(principalcomps)
    if pcs.ndim != 1 or not pcs.size or not np.all(np.isfinite(pcs)) or np.any(pcs != np.floor(pcs)) or np.any((pcs < 1) | (pcs > maps.shape[1])) or np.unique(pcs).size != pcs.size:
        raise ValueError("principalcomps must contain unique available one-based components")
    return maps, pcs.astype(int)


def coordinates(mni_coords: object, n_sources: int) -> np.ndarray:
    if mni_coords is None:
        path = files("broadness").joinpath("assets/MNI152_8mm_coord_dyi.mat")
        result = np.asarray(loadmat(path)["MNI8"], dtype=float)
    else:
        result = np.asarray(mni_coords, dtype=float)
    if result.shape != (n_sources, 3):
        raise ValueError("mni_coords must be sources × 3; provide custom coordinates for non-MNI data")
    return result


def prepare_maps(maps: np.ndarray, pcs: np.ndarray, coords: np.ndarray, normalization: str, distance: float | None) -> tuple[np.ndarray, np.ndarray, list[list[int]], float]:
    valid = np.all(np.isfinite(coords), axis=1)
    if valid.sum() < 2:
        raise ValueError("at least two valid MNI sources are required")
    analysis = maps[valid][:, pcs - 1, :, :].copy()
    if not np.all(np.isfinite(analysis)):
        raise ValueError("maps at valid MNI sources must be finite")
    if normalization == "rms":
        rms = np.sqrt(np.mean(analysis**2, axis=0, keepdims=True))
        if np.any(rms == 0):
            raise ValueError("cannot RMS-normalize an all-zero map")
        analysis /= rms
    elif normalization != "none":
        raise ValueError("mapnormalization must be 'none' or 'rms'")
    spatial = coords[valid]
    if distance is None:
        steps = np.concatenate([np.diff(np.unique(spatial[:, d])) for d in range(3)])
        steps = steps[steps > 1e-6]
        if not steps.size:
            raise ValueError("cannot infer MNI grid spacing")
        distance = float(steps.min())
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError("neighbourdistance must be positive")
    tree = cKDTree(spatial)
    neighbours = tree.query_ball_point(spatial, distance * 1.001)
    neighbours = [[j for j in row if j != i and np.sum((spatial[i] - spatial[j]) ** 2) > 1e-8] for i, row in enumerate(neighbours)]
    if not any(neighbours):
        raise ValueError("no spatial neighbours found; check coordinates or neighbourdistance")
    return valid, analysis, neighbours, distance


def label_clusters(stat: np.ndarray, threshold: np.ndarray | float, neighbours: list[list[int]], two_sided: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    stat = np.asarray(stat).ravel()
    threshold = np.broadcast_to(threshold, stat.shape)
    labels = np.zeros(stat.size, dtype=int)
    masses: list[float] = []
    signs: list[int] = []
    for sign in ((1, -1) if two_sided else (1,)):
        active = stat * sign >= threshold
        visited = np.zeros(stat.size, dtype=bool)
        for seed in np.flatnonzero(active):
            if visited[seed]:
                continue
            stack = [int(seed)]
            visited[seed] = True
            members = []
            while stack:
                current = stack.pop()
                members.append(current)
                for adj in neighbours[current]:
                    if active[adj] and not visited[adj]:
                        visited[adj] = True
                        stack.append(adj)
            labels[members] = len(masses) + 1
            masses.append(float(np.sum(np.abs(stat[members]))))
            signs.append(sign)
    return labels, np.asarray(masses), np.asarray(signs)


def cluster_outputs(stat: np.ndarray, threshold: np.ndarray | float, neighbours: list[list[int]], null: np.ndarray, valid: np.ndarray, coords: np.ndarray, alpha: float, two_sided: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Result]]:
    labels, masses, signs = label_clusters(stat, threshold, neighbours, two_sided)
    p_map = np.full(coords.shape[0], np.nan)
    significant = np.zeros(coords.shape[0], dtype=bool)
    full_labels = np.zeros(coords.shape[0], dtype=int)
    full_labels[valid] = labels
    source_ids = np.flatnonzero(valid)
    rows = []
    for i, (mass, sign) in enumerate(zip(masses, signs), start=1):
        p = (1 + np.count_nonzero(null >= mass)) / (null.size + 1)
        members = source_ids[labels == i]
        p_map[members] = p
        significant[members] = p <= alpha
        rows.append(Result(Cluster=i, Sign=int(sign), Mass=mass, CorrectedP=p,
                           Significant=bool(p <= alpha), SourceIndices=members + 1,
                           MNI_Coordinates=coords[members]))
    return full_labels, p_map, significant, rows


def paired_t(difference: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = difference.shape[-1]
    mean = difference.mean(axis=-1)
    sd = difference.std(axis=-1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = mean / (sd / np.sqrt(n))
        dz = mean / sd
    return np.where(np.isfinite(t), t, 0), mean, np.where(np.isfinite(dz), dz, np.nan)


def repeated_f(data: np.ndarray) -> np.ndarray:
    # sources × conditions × participants
    n_conditions, n_participants = data.shape[1:]
    grand = data.mean(axis=(1, 2), keepdims=True)
    condition_mean = data.mean(axis=2, keepdims=True)
    participant_mean = data.mean(axis=1, keepdims=True)
    condition_ss = n_participants * np.sum((condition_mean - grand) ** 2, axis=(1, 2))
    participant_ss = n_conditions * np.sum((participant_mean - grand) ** 2, axis=(1, 2))
    total_ss = np.sum((data - grand) ** 2, axis=(1, 2))
    error_ss = np.maximum(0, total_ss - condition_ss - participant_ss)
    with np.errstate(divide="ignore", invalid="ignore"):
        values = (condition_ss / (n_conditions - 1)) / (error_ss / ((n_participants - 1) * (n_conditions - 1)))
    return np.where(np.isfinite(values), values, 0)


def welch_t(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n1, n2 = x.shape[-1], y.shape[-1]
    m1, m2 = x.mean(axis=-1), y.mean(axis=-1)
    v1, v2 = x.var(axis=-1, ddof=1), y.var(axis=-1, ddof=1)
    a, b = v1 / n1, v2 / n2
    delta = m1 - m2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = delta / np.sqrt(a + b)
        df = (a + b) ** 2 / (a**2 / (n1 - 1) + b**2 / (n2 - 1))
        pooled = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
        g = (1 - 3 / (4 * (n1 + n2) - 9)) * delta / pooled
    return np.where(np.isfinite(t), t, 0), np.where(np.isfinite(df), df, 1), delta, np.where(np.isfinite(g), g, np.nan)


def welch_f(groups: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    n = np.asarray([group.shape[-1] for group in groups])
    means = np.stack([group.mean(axis=-1) for group in groups], axis=-1)
    variances = np.stack([group.var(axis=-1, ddof=1) for group in groups], axis=-1)
    invalid = np.any((variances <= 0) | ~np.isfinite(variances), axis=-1)
    weights = n / np.maximum(variances, np.finfo(float).tiny)
    weight_sum = weights.sum(axis=-1)
    weighted_mean = (weights * means).sum(axis=-1) / weight_sum
    correction = np.sum((1 - weights / weight_sum[..., None]) ** 2 / (n - 1), axis=-1)
    k = len(groups)
    df2 = (k**2 - 1) / (3 * correction)
    f = (np.sum(weights * (means - weighted_mean[..., None]) ** 2, axis=-1) / (k - 1)) / (1 + 2 * (k - 2) * correction / (k**2 - 1))
    return np.where(invalid | ~np.isfinite(f), 0, f), np.where(invalid | ~np.isfinite(df2), 1, df2)


SUMMARY_NAMES = ("MapStrength", "CentroidX", "CentroidY", "CentroidZ", "SpatialDispersion", "Lateralization")


def summaries(raw: np.ndarray, coords: np.ndarray) -> np.ndarray:
    # sources × selected components × conditions × participants
    weights = np.abs(raw)
    total = weights.sum(axis=0)
    strength = np.sqrt(np.mean(raw**2, axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        centroid = np.einsum("s...,sd->d...", weights, coords) / total
        second_moment = np.einsum("s...,s->...", weights, np.sum(coords**2, axis=1)) / total
        dispersion = np.sqrt(np.maximum(0, second_moment - np.sum(centroid**2, axis=0)))
        left = weights[coords[:, 0] < 0].sum(axis=0)
        right = weights[coords[:, 0] > 0].sum(axis=0)
        lateralization = (left - right) / (left + right)
    return np.concatenate((strength[None], centroid, dispersion[None], lateralization[None]), axis=0)


def fdr(values: np.ndarray, alpha: float) -> tuple[np.ndarray, np.ndarray, float]:
    return BROADNESS_FDRCorrection(values, alpha)


def validate_options(permutations: int, clusteralpha: float, alpha: float, randomseed: int, figuremode: str) -> None:
    if not isinstance(permutations, (int, np.integer)) or permutations < 3:
        raise ValueError("permutations must be an integer >= 3")
    if not 0 < alpha < 1 or not 0 < clusteralpha < 1:
        raise ValueError("alpha and clusteralpha must be between 0 and 1")
    if not isinstance(randomseed, (int, np.integer)):
        raise ValueError("randomseed must be an integer")
    if figuremode not in ("off", "show", "save", "both"):
        raise ValueError("figuremode must be 'off', 'show', 'save', or 'both'")


def pairs_or_default(pairs: object, n: int) -> np.ndarray:
    result = np.asarray(list(combinations(range(1, n + 1), 2)) if pairs is None else pairs)
    if result.ndim != 2 or result.shape[1] != 2 or not result.size or not np.all(np.isfinite(result)) or np.any(result != np.floor(result)) or np.any((result < 1) | (result > n)) or np.any(result[:, 0] == result[:, 1]):
        raise ValueError("invalid condition pairs")
    return result.astype(int)


def names_or_default(names: object, n: int, prefix: str) -> list[str]:
    if names is None:
        return [f"{prefix} {i}" for i in range(1, n + 1)]
    result = [str(x) for x in names]
    if len(result) != n:
        raise ValueError(f"provide exactly {n} {prefix.lower()} names")
    return result
