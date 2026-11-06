# Tutorial 2: compare conditions within participants

This follows MATLAB's `BROADNESS_Statistics_WithinParticipants_Example.m`.
Open [`02_BROADNESS_Statistics_WithinParticipants_Example.py`](../../examples/02_BROADNESS_Statistics_WithinParticipants_Example.py)
in VS Code after completing the [first-steps guide](../VS_CODE_FIRST_STEPS.md).

## Prepare the data

Download `DataReduced_AveragedOverParticipants_Example.mat` and at least four
`SUBJ_*.mat` files from the [Zenodo record](https://zenodo.org/records/17048137).
The participants are inside large ZIP archives. Put the extracted `.mat` files
directly in the repository's top-level `Data/` folder; the
[folder diagram](../VS_CODE_FIRST_STEPS.md#6-download-real-example-data-optional)
shows exactly where. The averaged file supplies the 776-sample time vector.

## Run from top to bottom

1. Run section **0**. Leave `TUTORIAL_MODE = True`, `MAX_PARTICIPANTS = 4`,
   `CORRECTION = "FDR"`, and `FIGURE_MODE = "show"` for the first run. If your
   files are elsewhere, change `data_path` in this section.
2. Run section **1**. It prints the participant file order and makes a single
   common PCA decomposition. In Variables, `network_time_series` has shape
   **time × networks × conditions × participants**; the initial lesson should
   show `(160, 2, 5, 4)`.
3. Run section **2**. It compares the same participants across conditions using
   paired tests. `CONDITION_PAIRS = ((1, 2),)` starts with one readable
   comparison; use `None` for all ten pairs. `CONDITION_STATS.PairedTests`
   contains p-values, FDR-adjusted p-values, significance, t statistics, and
   Cohen's *dz*. `CONDITION_STATS.ANOVA` tests the overall condition effect.
4. Run section **3**. The figures show condition averages. The `A` line refers
   to ANOVA; numbered lines refer to the condition pairs printed in the window.
5. Sections **4–6** are optional. Section 4 tests spatial maps and is off by
   default because it resamples across sources. In section 5, `BEHAVIOUR`
   must be one real score per participant in the **printed file order**.
   Section 6 correlates the same scores with RQA or phase-space measures.

If you change a setting in section 0 or 1, restart the Interactive Window and
run the cells in order again. `# %%` cells share the same Variables explorer
only while that kernel remains open.

## Interpreting this lesson

Tutorial mode uses 128 sources and 160 time samples from the real files, with
four participants. It demonstrates the workflow, **not** a defensible
scientific result. `CORRECTION = "cluster"` calls the temporal permutation
functions; the 100 tutorial permutations have coarse p-value resolution.
Full inference needs a suitable sample, the full resolution
(`TUTORIAL_MODE = False`), an intentional participant selection, and enough
permutations. Set `MAX_PARTICIPANTS = None` only when you are ready to load all
available files. Spatial statistics can take much longer than the first three
sections. Corrected cluster p-values apply to whole clusters, not to each
individual time sample or source. The loader stops before allocating a data
array larger than 2 GiB unless you explicitly set `ALLOW_LARGE = True` after
checking available memory.
The loader also stops before reading more than 20 participant files by default,
even in tutorial mode, because those files require substantial disk I/O.
