# Python porting status

Reference MATLAB version: `main` at `7a8d54b` (2 November 2026), now copied to
the sibling `matlab` directory. Compare against that version before implementing
or updating each Python function.

The Python package currently has 12 of the 20 functions in
`matlab/BROADNESS_Functions`, plus `BROADNESS_Startup`. “Implemented” means
source and Python tests are present. Direct MATLAB numerical comparisons are
still pending; MATLAB's command-line launch did not complete in this workspace.

| MATLAB function | Python status | Main work remaining |
| --- | --- | --- |
| `BROADNESS_AlternativeNetworkEstimation_ICA` | Implemented | JADE decomposition and MATLAB output fields; direct MATLAB numerical comparison pending. |
| `BROADNESS_BehaviourCorrelation` | Implemented | Pearson/Spearman correlations, pairwise valid participants, and FDR. |
| `BROADNESS_EffectiveDimensionality` | Existing | Check against the current MATLAB behavior. |
| `BROADNESS_FDRCorrection` | Implemented | Benjamini–Hochberg correction and NaN handling. |
| `BROADNESS_FigureSettings` | Missing | Common plot and export options. |
| `BROADNESS_FinalizeFigure` | Missing | Common figure export. |
| `BROADNESS_NetworkEstimation` | Updated | `components_num`, common PCA centering, retained components, and participant spatial patterns have Python tests; direct MATLAB comparison pending. |
| `BROADNESS_PhaseSpaceStatistics` | Missing | Phase-space tests, metrics, and corrections. |
| `BROADNESS_PhaseSpace_RQA` | Missing | Recurrence metrics and plots. |
| `BROADNESS_Plot_ActivationTimeseries` | Missing | Statistical time-series plots. |
| `BROADNESS_RepeatedMeasuresANOVA` | Implemented | Time-resolved condition test, Greenhouse–Geisser correction, and family-wide FDR. |
| `BROADNESS_SignificantSamplesToWindows` | Implemented | Convert masks into time windows. |
| `BROADNESS_SpatialActivationClustering` | Missing | Clustering, results, plots, and optional NIfTI export. |
| `BROADNESS_SpatialPatternGroupStatistics` | Missing | Independent-group spatial inference. |
| `BROADNESS_SpatialPatternStatistics` | Missing | Within-participant spatial inference. |
| `BROADNESS_TemporalClusterCorrection` | Implemented | Paired whole-series sign flips; cluster-level inference per network and condition pair. |
| `BROADNESS_TemporalClusterCorrection_ANOVA` | Implemented | Complete condition-label permutations within participants; per-network nulls. |
| `BROADNESS_TemporalClusterCorrection_GroupANOVA` | Implemented | Welch ANOVA with group-centered wild bootstrap; approximate cluster inference. |
| `BROADNESS_TemporalClusterCorrection_Groups` | Implemented | Pairwise Welch tests with group-centered wild bootstrap; approximate cluster inference. |
| `BROADNESS_Visualizer` | Existing, needs update | New figure mode, layout, export options, and output structure. |

The three requested examples are all pending: `BROADNESS_Example`,
`BROADNESS_Statistics_WithinParticipants_Example`, and
`BROADNESS_Statistics_MultipleGroupsParticipants_Example`.

## External functionality

| MATLAB external component | Python approach |
| --- | --- |
| `jadeR.m` | Python now implements real-valued JADE with PCA whitening, fourth-order cumulants, joint diagonalization, and MATLAB's ordering/sign conventions. Direct MATLAB comparison remains pending. |
| `dl.m`, `tt.m`, `entropy.m` | Implement the recurrence line statistics and entropy used by RQA. |
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
