# Tutorial 3: compare independent participant groups

This follows MATLAB's `BROADNESS_Statistics_MultipleGroupsParticipants_Example.m`.
Open [`03_BROADNESS_Statistics_MultipleGroupsParticipants_Example.py`](../../examples/03_BROADNESS_Statistics_MultipleGroupsParticipants_Example.py)
after completing the [VS Code first steps](../VS_CODE_FIRST_STEPS.md).

## Prepare the data

From the [Zenodo record](https://zenodo.org/records/17048137), obtain the
averaged `.mat` file and at least six participant `SUBJ_*.mat` files. Put the
extracted files directly in the repository's top-level `Data/` folder.
Participant files come from large ZIP archives. The averaged file supplies the
time vector. See the [exact folder layout](../VS_CODE_FIRST_STEPS.md#6-download-real-example-data-optional).

**The participant `.mat` files contain source data, not the group assignment
needed for this example.** In tutorial mode, the first half of the selected files are assigned
to “Demonstration A” and the rest to “Demonstration B.” Those labels are
arbitrary. Their p-values and effect sizes must not be read as a scientific
comparison.

## Run from top to bottom

1. Run section **0**. Keep `TUTORIAL_MODE = True`, `MAX_PARTICIPANTS = 6`,
   `CORRECTION = "FDR"`, and `FIGURE_MODE = "show"` for your first run.
2. Run section **1**. Read the printed participant order. The initial lesson
   should produce `network_time_series` with shape `(160, 2, 5, 6)`.
3. Run section **2**. It estimates one common decomposition and compares
   independent groups using Welch tests. Inspect
   `GROUP_STATS.BetweenGroupComparisons` for p-values, significance, t
   statistics, and Hedges' *g*. With three or more groups it also runs an
   overall Welch ANOVA.
4. Run section **3**. Inspect the printed significance-line key and the group
   time-series plots.
5. Sections **4–6** are optional. Section 4 fits a condition × group mixed
   model in a time window chosen in advance; it requires **real group labels**
   and is off initially. Section 5 adds a behavioural score to that model;
   section 6 tests participant-specific spatial maps. Section 6 can run with
   demonstration groups to check software behavior, but its results remain
   arbitrary.

## Replace demonstration labels with real groups

In section 0, enter `GROUP_LABELS_BY_FILE` with one entry for **every** selected
participant file, for example:

```python
GROUP_LABELS_BY_FILE = {
    "SUBJ_0002.mat": 1,
    "SUBJ_0007.mat": 2,
    # Continue for every file printed by section 1.
}
GROUP_NAMES = ["Your group 1", "Your group 2"]
```

The mapping above only illustrates the format; it does **not** assert the true
group membership of those subjects. The script checks for missing assignments
and requires at least three participants in each group. Set
`TUTORIAL_MODE = False` to use all sources and time samples; set
`MAX_PARTICIPANTS = None` only when you intend to load every participant file
in `Data/`. Statistical inference also requires a suitable study design and
sample size. Group cluster correction uses an approximate wild bootstrap;
the 100 tutorial resamples are only for checking software behavior. The
Python mixed-model fit is a close analogue of the MATLAB example, but its
numerical results and inferential method have not yet been compared directly
against MATLAB. A mixed model can be singular with small or low-variance groups;
the script reports that fit failure rather than treating it as a result. The
loader stops before allocating a data array larger than
2 GiB unless you explicitly set `ALLOW_LARGE = True` after checking available
memory. It also stops before reading more than 20 participant files by default,
even in tutorial mode, because those files require substantial disk I/O.
