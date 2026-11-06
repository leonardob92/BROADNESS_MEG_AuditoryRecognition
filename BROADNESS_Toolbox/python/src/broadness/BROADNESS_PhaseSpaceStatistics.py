"""Descriptive and inferential statistics for BROADNESS phase-space paths."""

from __future__ import annotations

from itertools import combinations
from types import SimpleNamespace as Result

import numpy as np
from numpy.typing import ArrayLike
from scipy import stats

from .BROADNESS_FDRCorrection import BROADNESS_FDRCorrection
from ._spatial import field


METRIC_NAMES = ("Dispersion", "Speed", "Acceleration", "BaselineDisplacement", "CumulativePathLength", "TotalPathLength")


def _nanmean(values: np.ndarray, axis: int) -> np.ndarray:
    valid = np.isfinite(values)
    total = np.sum(np.where(valid, values, 0), axis=axis)
    count = np.sum(valid, axis=axis)
    return np.divide(total, count, out=np.full_like(total, np.nan, dtype=float), where=count > 0)


def _adjust(values: np.ndarray, alpha: float, fdr: bool) -> np.ndarray:
    if fdr:
        return BROADNESS_FDRCorrection(values, alpha)[1]
    return values.copy()


def _groups(groupindices: object, groupnames: object, n: int) -> tuple[list[np.ndarray], list[str], np.ndarray]:
    if groupindices is None:
        members = [np.arange(n)]
    elif isinstance(groupindices, (list, tuple)) and len(groupindices) and not np.isscalar(groupindices[0]):
        members = [np.asarray(row, dtype=int).ravel() - 1 for row in groupindices]
    else:
        labels = np.asarray(groupindices, dtype=float).ravel()
        if labels.size != n:
            raise ValueError("groupindices must be group member lists or one label per participant")
        levels = list(dict.fromkeys(labels[np.isfinite(labels)]))
        members = [np.flatnonzero(labels == level) for level in levels]
    if not members or any(len(row) == 0 or np.any((row < 0) | (row >= n)) for row in members):
        raise ValueError("groupindices contain invalid participants")
    if np.unique(np.concatenate(members)).size != sum(map(len, members)):
        raise ValueError("participants cannot occur in more than one group")
    names = (["AllParticipants"] if len(members) == 1 else [f"Group{i}" for i in range(1, len(members) + 1)]) if groupnames is None else [str(x) for x in groupnames]
    if len(names) != len(members):
        raise ValueError("groupnames must match groupindices")
    labels = np.full(n, np.nan)
    for index, row in enumerate(members, start=1):
        labels[row] = index
    return members, names, labels


def _hotelling(a: np.ndarray, b: np.ndarray, paired: bool) -> Result:
    # time × dimensions × participants
    n_time = a.shape[0]
    output = Result(**{key: np.full(n_time, np.nan) for key in ("T2", "F", "DF1", "DF2", "P", "EffectiveDimensions")})
    for t in range(n_time):
        x, y = a[t].T, b[t].T
        if paired:
            difference = x - y
            difference = difference[np.all(np.isfinite(difference), axis=1)]
            n = len(difference)
            if n < 2:
                continue
            covariance = np.atleast_2d(np.cov(difference, rowvar=False))
            mean = difference.mean(axis=0)
            factor = n
            denominator = n - 1
        else:
            x = x[np.all(np.isfinite(x), axis=1)]
            y = y[np.all(np.isfinite(y), axis=1)]
            n1, n2 = len(x), len(y)
            if min(n1, n2) < 2:
                continue
            covariance = ((n1 - 1) * np.atleast_2d(np.cov(x, rowvar=False)) + (n2 - 1) * np.atleast_2d(np.cov(y, rowvar=False))) / (n1 + n2 - 2)
            mean = x.mean(axis=0) - y.mean(axis=0)
            factor = n1 * n2 / (n1 + n2)
            denominator = n1 + n2 - 2
        if not np.all(np.isfinite(covariance)):
            continue
        rank = np.linalg.matrix_rank(covariance)
        df2 = (n - rank) if paired else (n1 + n2 - rank - 1)
        if rank < 1 or df2 < 1:
            continue
        t2 = float(factor * mean @ np.linalg.pinv(covariance) @ mean)
        f = df2 * t2 / (rank * denominator)
        for key, value in (("T2", t2), ("F", f), ("DF1", rank), ("DF2", df2),
                           ("P", stats.f.sf(f, rank, df2)), ("EffectiveDimensions", rank)):
            getattr(output, key)[t] = value
    return output


def _hotelling_comparison(a: np.ndarray, b: np.ndarray, time: np.ndarray, level1: str, level2: str, context: str, kind: str, paired: bool) -> Result:
    mean_a = _nanmean(a, axis=2)
    mean_b = _nanmean(b, axis=2)
    difference = mean_a - mean_b
    return Result(Label=f"{level1}_vs_{level2}__{context}", Comparison=kind,
                  Level1=level1, Level2=level2, Context=context, Time=time,
                  Level1MeanCoordinates=mean_a, Level2MeanCoordinates=mean_b,
                  SignedComponentDifference=difference,
                  EuclideanDistance=np.linalg.norm(difference, axis=1),
                  Hotelling=_hotelling(a, b, paired))


def _scalar_comparison(a: np.ndarray, b: np.ndarray, time: np.ndarray, level1: str, level2: str, context: str, kind: str, paired: bool, testtype: str) -> Result:
    comparison = Result(Label=f"{level1}_vs_{level2}__{context}", Comparison=kind,
                        Level1=level1, Level2=level2, Context=context, Time=time,
                        Level1Mean=_nanmean(a, axis=1), Level2Mean=_nanmean(b, axis=1))
    comparison.MeanDifference = comparison.Level1Mean - comparison.Level2Mean
    for key in (("Parametric", "Nonparametric") if testtype == "both" else (("Parametric",) if testtype == "parametric" else ("Nonparametric",))):
        output = Result(Statistic=np.full(len(time), np.nan), DF=np.full(len(time), np.nan),
                        P=np.full(len(time), np.nan), P_FDR=np.full(len(time), np.nan),
                        Significant=np.zeros(len(time), dtype=bool),
                        N1=np.full(len(time), np.nan), N2=np.full(len(time), np.nan))
        for t in range(len(time)):
            x, y = a[t], b[t]
            if paired:
                valid = np.isfinite(x) & np.isfinite(y)
                x, y = x[valid], y[valid]
            else:
                x, y = x[np.isfinite(x)], y[np.isfinite(y)]
            output.N1[t], output.N2[t] = len(x), len(y)
            if min(len(x), len(y)) < 2:
                continue
            if key == "Parametric":
                result = stats.ttest_rel(x, y) if paired else stats.ttest_ind(x, y, equal_var=False)
                output.Statistic[t], output.P[t] = result.statistic, result.pvalue
                output.DF[t] = result.df
            elif paired:
                if np.all(x == y):
                    output.Statistic[t], output.P[t] = 0, 1
                else:
                    result = stats.wilcoxon(x, y)
                    output.Statistic[t], output.P[t] = result.statistic, result.pvalue
            else:
                result = stats.mannwhitneyu(x, y, alternative="two-sided")
                output.Statistic[t], output.P[t] = result.statistic, result.pvalue
        setattr(comparison, key, output)
    return comparison


def _family(comparisons: list[Result], key: str, alpha: float, use_fdr: bool) -> None:
    if not comparisons:
        return
    p = np.concatenate([getattr(row, key).P for row in comparisons])
    adjusted = _adjust(p, alpha, use_fdr)
    start = 0
    for row in comparisons:
        result = getattr(row, key)
        count = len(result.P)
        result.P_FDR = adjusted[start:start + count]
        result.Significant = result.P_FDR < alpha
        start += count


def _mixed_anova(values: np.ndarray, labels: np.ndarray, time: np.ndarray, alpha: float, use_fdr: bool) -> Result:
    n_time, n_conditions, _ = values.shape
    effects = {name: Result(**{key: np.full(n_time, np.nan) for key in ("F", "DF1", "DF2", "P", "EpsilonGG", "DF1_GG", "DF2_GG", "P_GG", "P_FDR")}, Significant=np.zeros(n_time, dtype=bool)) for name in ("Condition", "Group", "Interaction")}
    for t in range(n_time):
        cell = values[t].T
        valid = np.all(np.isfinite(cell), axis=1) & np.isfinite(labels)
        cell, groups = cell[valid], labels[valid]
        n = len(cell)
        if n < 2:
            continue
        levels = np.unique(groups)
        k = len(levels)
        grand = cell.mean()
        participant_mean = cell.mean(axis=1)
        condition_mean = cell.mean(axis=0)
        condition_ss = n * np.sum((condition_mean - grand) ** 2)
        if k == 1:
            if n_conditions < 2:
                continue
            error = cell - participant_mean[:, None] - condition_mean + grand
            ss = {"Condition": (condition_ss, n_conditions - 1, np.sum(error**2), (n - 1) * (n_conditions - 1))}
            residual = cell - cell.mean(axis=0)
        else:
            group_means = np.array([cell[groups == level].mean() for level in levels])
            group_condition = np.array([cell[groups == level].mean(axis=0) for level in levels])
            sizes = np.array([np.count_nonzero(groups == level) for level in levels])
            subject_ss = n_conditions * sum(np.sum((participant_mean[groups == level] - group_means[i]) ** 2) for i, level in enumerate(levels))
            group_ss = n_conditions * np.sum(sizes * (group_means - grand) ** 2)
            interaction_ss = sum(sizes[i] * np.sum((group_condition[i] - group_means[i] - condition_mean + grand) ** 2) for i in range(k))
            error_ss = sum(np.sum((cell[groups == level] - participant_mean[groups == level, None] - group_condition[i] + group_means[i]) ** 2) for i, level in enumerate(levels))
            ss = {"Group": (group_ss, k - 1, subject_ss, n - k),
                  "Condition": (condition_ss, n_conditions - 1, error_ss, (n - k) * (n_conditions - 1)),
                  "Interaction": (interaction_ss, (k - 1) * (n_conditions - 1), error_ss, (n - k) * (n_conditions - 1))}
            residual = cell.copy()
            for level in levels:
                residual[groups == level] -= cell[groups == level].mean(axis=0)
        epsilon = 1.0
        if n_conditions > 2 and n > 1:
            covariance = np.cov(residual, rowvar=False)
            centering = np.eye(n_conditions) - np.ones((n_conditions, n_conditions)) / n_conditions
            centered = centering @ covariance @ centering
            denom = (n_conditions - 1) * np.trace(centered @ centered)
            if denom > 0:
                epsilon = float(np.clip(np.trace(centered) ** 2 / denom, 1 / (n_conditions - 1), 1))
        for name, (effect_ss, df1, error_ss, df2) in ss.items():
            if df1 < 1 or df2 < 1 or error_ss < 0:
                continue
            f = ((effect_ss / df1) / (error_ss / df2)) if error_ss else (np.inf if effect_ss else 0.0)
            output = effects[name]
            for key, value in (("F", f), ("DF1", df1), ("DF2", df2), ("P", stats.f.sf(f, df1, df2))):
                getattr(output, key)[t] = value
            correction = epsilon if name != "Group" else 1.0
            output.EpsilonGG[t] = correction
            output.DF1_GG[t], output.DF2_GG[t] = df1 * correction, df2 * correction
            output.P_GG[t] = stats.f.sf(f, df1 * correction, df2 * correction)
    for output in effects.values():
        output.P_FDR = _adjust(output.P_GG, alpha, use_fdr)
        output.Significant = output.P_FDR < alpha
    return Result(Time=time, **effects)


def BROADNESS_PhaseSpaceStatistics(
    RQA_BROADNESS: object, *, conditionindices: ArrayLike | None = None,
    conditionnames: list[str] | None = None, groupindices: object = None,
    groupnames: list[str] | None = None, statistics: str = "on",
    testtype: str = "parametric", pairwise: str = "on", omnibus: str = "on",
    alpha: float = 0.05, fdr: str = "on", baseline: ArrayLike | None = None,
    smoothingsamples: int = 1, testmetrics: list[str] | tuple[str, ...] = METRIC_NAMES,
) -> Result:
    """Calculate trajectory metrics, Hotelling tests, and scalar tests.

    ``groupindices`` accepts one-based participant lists or one group label
    per participant. Inference uses complete observations at each time point.
    """
    phase = field(RQA_BROADNESS, "PhaseSpace")
    time = np.asarray(field(phase, "Time"), dtype=float).ravel()
    pcs = np.asarray(field(phase, "PCs"), dtype=int).ravel()
    cells = np.asarray(field(phase, "ParticipantCoordinates"), dtype=object)
    if cells.ndim != 2 or not cells.size or time.size < 3 or not np.all(np.isfinite(time)) or not np.all(np.diff(time) > 0):
        raise ValueError("PhaseSpace requires increasing Time and a conditions × participants coordinate array")
    n_original, n_participants = cells.shape
    selected = np.arange(1, n_original + 1) if conditionindices is None else np.asarray(conditionindices)
    if selected.ndim != 1 or not selected.size or np.any((selected < 1) | (selected > n_original)) or np.unique(selected).size != selected.size or np.any(selected != np.floor(selected)):
        raise ValueError("conditionindices must contain unique valid one-based conditions")
    selected = selected.astype(int)
    names = [f"Condition{i}" for i in range(1, len(selected) + 1)] if conditionnames is None else [str(x) for x in conditionnames]
    if len(names) != len(selected):
        raise ValueError("conditionnames must match selected conditions")
    members, group_names, labels = _groups(groupindices, groupnames, n_participants)
    included = np.unique(np.concatenate(members))
    first = np.asarray(cells[selected[0] - 1, 0], dtype=float)
    if first.ndim != 2 or first.shape != (len(time), len(pcs)):
        raise ValueError("trajectory dimensions must be time × selected components")
    trajectories = np.empty((len(time), len(pcs), len(selected), n_participants))
    for condition, old_index in enumerate(selected):
        for participant in range(n_participants):
            value = np.asarray(cells[old_index - 1, participant], dtype=float)
            if value.shape != first.shape:
                raise ValueError("all trajectories must have the same shape")
            trajectories[:, :, condition, participant] = value
    if not isinstance(smoothingsamples, int) or smoothingsamples < 1:
        raise ValueError("smoothingsamples must be a positive integer")
    if smoothingsamples > 1:
        previous = trajectories.copy()
        left = (smoothingsamples - 1) // 2
        right = smoothingsamples - left - 1
        for t in range(len(time)):
            trajectories[t] = np.mean(previous[max(0, t - left):min(len(time), t + right + 1)], axis=0)
    if baseline is None:
        baseline_mask = time <= 0
        if not np.any(baseline_mask):
            baseline_mask[0] = True
    else:
        bounds = np.asarray(baseline, dtype=float).ravel()
        if bounds.size != 2 or not np.all(np.isfinite(bounds)):
            raise ValueError("baseline must contain two finite times")
        lo, hi = sorted(bounds)
        baseline_mask = (time >= lo) & (time <= hi)
        if not np.any(baseline_mask):
            raise ValueError("baseline contains no time samples")
    for option, value in (("statistics", statistics), ("pairwise", pairwise), ("omnibus", omnibus), ("fdr", fdr)):
        if value not in ("on", "off"):
            raise ValueError(f"{option} must be 'on' or 'off'")
    if testtype not in ("parametric", "nonparametric", "both") or not 0 < alpha < 1:
        raise ValueError("invalid testtype or alpha")
    metric_names = list(dict.fromkeys(next((valid for valid in METRIC_NAMES if valid.lower() == str(name).lower()), None) for name in ([testmetrics] if isinstance(testmetrics, str) else testmetrics)))
    if None in metric_names:
        raise ValueError("testmetrics contains an unknown metric")
    n_time, n_dim, n_conditions, _ = trajectories.shape
    mean_coordinates = _nanmean(trajectories[:, :, :, included], axis=3)
    velocity = np.full_like(trajectories, np.nan)
    acceleration_vector = np.full_like(trajectories, np.nan)
    velocity[1:] = np.diff(trajectories, axis=0) / np.diff(time)[:, None, None, None]
    acceleration_vector[2:] = np.diff(velocity[1:], axis=0) / ((np.diff(time)[:-1] + np.diff(time)[1:]) / 2)[:, None, None, None]
    speed = np.linalg.norm(velocity, axis=1)
    acceleration = np.linalg.norm(acceleration_vector, axis=1)
    steps = np.linalg.norm(np.diff(trajectories, axis=0), axis=1)
    cumulative = np.concatenate((np.zeros((1, n_conditions, n_participants)), np.cumsum(steps, axis=0)))
    total = np.nansum(steps, axis=0)
    displacement = np.empty((n_time, n_conditions, n_participants))
    dispersion = np.empty_like(displacement)
    for participant in range(n_participants):
        others = included[included != participant] if participant in included else included
        centroid = _nanmean(_nanmean(trajectories[:, :, :, others], axis=3), axis=2) if len(others) else _nanmean(trajectories[:, :, :, participant], axis=2)
        for condition in range(n_conditions):
            path = trajectories[:, :, condition, participant]
            baseline_centroid = _nanmean(path[baseline_mask], axis=0)
            displacement[:, condition, participant] = np.linalg.norm(path - baseline_centroid, axis=1)
            dispersion[:, condition, participant] = np.linalg.norm(path - centroid, axis=1)
    comparisons = []
    for first_condition, second_condition in combinations(range(n_conditions), 2):
        difference = mean_coordinates[:, :, first_condition] - mean_coordinates[:, :, second_condition]
        comparisons.append(Result(Label=f"{names[first_condition]}_vs_{names[second_condition]}",
                                  Condition1=names[first_condition], Condition2=names[second_condition],
                                  Time=time, SignedComponentDifference=difference,
                                  EuclideanDistance=np.linalg.norm(difference, axis=1)))
    metrics = Result(Velocity=velocity, Speed=speed, AccelerationVector=acceleration_vector,
                     Acceleration=acceleration, Dispersion=dispersion,
                     BaselineDisplacement=displacement, CumulativePathLength=cumulative,
                     TotalPathLength=total)
    settings = Result(ConditionIndices=selected, ConditionNames=names,
                      GroupMembers=[row + 1 for row in members], GroupNames=group_names,
                      GroupLabels=labels, Baseline=np.array([time[baseline_mask][0], time[baseline_mask][-1]]),
                      SmoothingSamples=smoothingsamples, Alpha=alpha, FDR=fdr,
                      FDRFamily="Each comparison category across its time points and pairwise comparisons",
                      TestType=testtype)
    result = Result(Time=time, PCs=pcs, Trajectories=trajectories, Metrics=metrics,
                    Descriptive=Result(ConditionMeanCoordinates=mean_coordinates,
                                       CentroidComparisons=comparisons), Settings=settings)
    if statistics == "off" or len(included) < 2:
        result.Tests = Result(Enabled=False, Reason="Statistics disabled or fewer than two participants")
        return result
    test = Result(Enabled=True)
    contexts = list(zip(members, group_names))
    if len(members) > 1:
        contexts = [(included, "AllParticipants"), *contexts]
    categories = ("ConditionComparisons", "GroupComparisons", "Interactions")
    phase_tests = Result(**{name: [] for name in categories})
    if pairwise == "on":
        for ids, context in contexts:
            for a, b in combinations(range(n_conditions), 2):
                first = trajectories[:, :, a, ids]
                second = trajectories[:, :, b, ids]
                phase_tests.ConditionComparisons.append(_hotelling_comparison(first, second, time, names[a], names[b], context, "conditions", True))
        for condition in range(n_conditions):
            for a, b in combinations(range(len(members)), 2):
                phase_tests.GroupComparisons.append(_hotelling_comparison(
                    trajectories[:, :, condition, members[a]], trajectories[:, :, condition, members[b]],
                    time, group_names[a], group_names[b], names[condition], "groups", False))
        for a, b in combinations(range(len(members)), 2):
            for first, second in combinations(range(n_conditions), 2):
                difference_a = trajectories[:, :, first, members[a]] - trajectories[:, :, second, members[a]]
                difference_b = trajectories[:, :, first, members[b]] - trajectories[:, :, second, members[b]]
                context = f"{names[first]}_minus_{names[second]}"
                phase_tests.Interactions.append(_hotelling_comparison(difference_a, difference_b, time, group_names[a], group_names[b], context, "interaction", False))
        for category in categories:
            _family(getattr(phase_tests, category), "Hotelling", alpha, fdr == "on")
    test.PhaseSpace = phase_tests
    metric_tests = Result()
    for name in metric_names:
        data = total[None] if name == "TotalPathLength" else getattr(metrics, name)
        metric_time = time[-1:] if name == "TotalPathLength" else time
        output = Result(**{category: [] for category in categories}, Omnibus=Result())
        if pairwise == "on":
            for ids, context in contexts:
                for a, b in combinations(range(n_conditions), 2):
                    output.ConditionComparisons.append(_scalar_comparison(
                        data[:, a, ids], data[:, b, ids], metric_time,
                        names[a], names[b], context, "conditions", True, testtype))
            for condition in range(n_conditions):
                for a, b in combinations(range(len(members)), 2):
                    output.GroupComparisons.append(_scalar_comparison(
                        data[:, condition, members[a]], data[:, condition, members[b]], metric_time,
                        group_names[a], group_names[b], names[condition], "groups", False, testtype))
            for a, b in combinations(range(len(members)), 2):
                for first, second in combinations(range(n_conditions), 2):
                    x = data[:, first, members[a]] - data[:, second, members[a]]
                    y = data[:, first, members[b]] - data[:, second, members[b]]
                    output.Interactions.append(_scalar_comparison(x, y, metric_time,
                        group_names[a], group_names[b], f"{names[first]}_minus_{names[second]}",
                        "interaction", False, testtype))
            for category in categories:
                rows = getattr(output, category)
                for key in (("Parametric", "Nonparametric") if testtype == "both" else (("Parametric",) if testtype == "parametric" else ("Nonparametric",))):
                    _family(rows, key, alpha, fdr == "on")
        if omnibus == "on":
            output.Omnibus.Parametric = _mixed_anova(data, labels, metric_time, alpha, fdr == "on")
        setattr(metric_tests, name, output)
    test.Metrics = metric_tests
    result.Tests = test
    return result


__all__ = ["BROADNESS_PhaseSpaceStatistics"]
