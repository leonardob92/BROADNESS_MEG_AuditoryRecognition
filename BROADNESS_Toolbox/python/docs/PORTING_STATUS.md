# Python porting status

Reference MATLAB version: `main` at `7a8d54b` (2 November 2026), now copied to
the sibling `matlab` directory. Compare against that version before implementing
or updating each Python function.

The Python package now has corresponding modules for all 20 functions in
`matlab/BROADNESS_Functions`, plus `BROADNESS_Startup`. “Implemented” means
source and Python tests are present. Direct MATLAB numerical comparisons are
still pending; MATLAB's command-line launch did not complete in this workspace.

| MATLAB function | Python status | Main work remaining |
| --- | --- | --- |
| `BROADNESS_AlternativeNetworkEstimation_ICA` | Implemented | JADE decomposition and MATLAB output fields; direct MATLAB numerical comparison pending. |
| `BROADNESS_BehaviourCorrelation` | Implemented | Pearson/Spearman correlations, pairwise valid participants, and FDR. |
| `BROADNESS_EffectiveDimensionality` | Existing | Check against the current MATLAB behavior. |
| `BROADNESS_FDRCorrection` | Implemented | Benjamini–Hochberg correction and NaN handling. |
| `BROADNESS_FigureSettings` | Implemented | Shared mode, layout, format, prefix, and output-folder validation. |
| `BROADNESS_FinalizeFigure` | Implemented | PNG/PDF export and Python-native editable `.fig` pickle. Python `.fig` files are not MATLAB-compatible. |
| `BROADNESS_NetworkEstimation` | Updated | `components_num`, common PCA centering, retained components, and participant spatial patterns have Python tests; direct MATLAB comparison pending. |
| `BROADNESS_PhaseSpaceStatistics` | Implemented | Descriptive metrics, Hotelling tests, scalar tests, mixed ANOVA, and corrections; direct MATLAB numerical comparison pending. |
| `BROADNESS_PhaseSpace_RQA` | Implemented | Recurrence matrices, eight measures, individual/summary plots, and optional animation; direct MATLAB numerical comparison pending. |
| `BROADNESS_Plot_ActivationTimeseries` | Implemented | Condition/group means, standard error styles, significance windows, and figure layouts. |
| `BROADNESS_RepeatedMeasuresANOVA` | Implemented | Time-resolved condition test, Greenhouse–Geisser correction, and family-wide FDR. |
| `BROADNESS_SignificantSamplesToWindows` | Implemented | Convert masks into time windows. |
| `BROADNESS_SpatialActivationClustering` | Implemented | Thresholding, breadth ordering, silhouette selection, tables, figures, and optional NIfTI export. Python defaults use 10 k-means restarts and a 500-source silhouette sample to control cost. |
| `BROADNESS_SpatialPatternGroupStatistics` | Implemented | Welch maps, wild-bootstrap cluster inference, spatial summaries, FDR, and figures. |
| `BROADNESS_SpatialPatternStatistics` | Implemented | Paired sign-flip maps, optional repeated-measures F map, spatial summaries, FDR, and figures. |
| `BROADNESS_TemporalClusterCorrection` | Implemented | Paired whole-series sign flips; cluster-level inference per network and condition pair. |
| `BROADNESS_TemporalClusterCorrection_ANOVA` | Implemented | Complete condition-label permutations within participants; per-network nulls. |
| `BROADNESS_TemporalClusterCorrection_GroupANOVA` | Implemented | Welch ANOVA with group-centered wild bootstrap; approximate cluster inference. |
| `BROADNESS_TemporalClusterCorrection_Groups` | Implemented | Pairwise Welch tests with group-centered wild bootstrap; approximate cluster inference. |
| `BROADNESS_Visualizer` | Updated | MATLAB-style figure mode, layout, export options, and file listing; NIfTI export now uses `OutputPath/BROADNESS_nifti`. |

The three MATLAB-style real-data Python examples are now available in
[`examples/`](../examples/) with individual guides in
[`docs/tutorials/`](tutorials/). The general and two statistics scripts mirror
the main MATLAB section order, but tutorial mode samples sources and time to
limit first-run cost. The multiple-group tutorial uses arbitrary demonstration
labels until real group metadata is supplied. Direct MATLAB numerical parity
checks and full-data scientific validation remain pending.

## External functionality

| MATLAB external component | Python approach |
| --- | --- |
| `jadeR.m` | Python now implements real-valued JADE with PCA whitening, fourth-order cumulants, joint diagonalization, and MATLAB's ordering/sign conventions. Direct MATLAB comparison remains pending. |
| `dl.m`, `tt.m`, `entropy.m` | Python now computes diagonal/vertical recurrence runs and natural-log entropy directly for RQA. Direct MATLAB numerical comparison remains pending. |
| `NIfTI_20140122` | Use `nibabel` for the needed NIfTI reads and writes; the 42 MATLAB helper files do not need individual translations. |
| MNI coordinates and brain templates | The small 8-mm coordinate/template files and compact 1-mm surface mesh are packaged with Python. Verify the installed package can use them outside the repository. |

## Shared conventions and verification

- Preserve the MATLAB axis order: sources × time × conditions × participants
  for input data, and time × components × conditions × participants for network
  time series.
- Preserve one-based component and participant selections at the public API;
  internal NumPy indexing remains zero-based.
- Use small, deterministic fixtures for numerical comparisons and keep
  expensive permutation/bootstrap runs out of the default test suite.
- Provide ordinary `pip` installation and optional dependencies for plots,
  NIfTI files, and MATLAB data import where possible.
- Keep participant data and generated output out of Git.
