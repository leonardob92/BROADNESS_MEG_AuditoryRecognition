# Python porting status

Reference MATLAB version: `main` at `7a8d54b` (2 November 2026), now copied to
the sibling `matlab` directory. Compare against that version before implementing
or updating each Python function.

The Python package currently has three of the 20 functions in
`matlab/BROADNESS_Functions`, plus `BROADNESS_Startup`. “Existing” means that source and
tests are present; it does not mean parity with the newer MATLAB version.

| MATLAB function | Python status | Main work remaining |
| --- | --- | --- |
| `BROADNESS_AlternativeNetworkEstimation_ICA` | Missing | ICA decomposition and output contract; assess JADE compatibility. |
| `BROADNESS_BehaviourCorrelation` | Missing | Pearson/Spearman correlations and FDR. |
| `BROADNESS_EffectiveDimensionality` | Existing | Check against the current MATLAB behavior. |
| `BROADNESS_FDRCorrection` | Missing | Benjamini–Hochberg correction and NaN handling. |
| `BROADNESS_FigureSettings` | Missing | Common plot and export options. |
| `BROADNESS_FinalizeFigure` | Missing | Common figure export. |
| `BROADNESS_NetworkEstimation` | Existing, needs update | `components_num`, common PCA centering, retained components, and participant spatial patterns. |
| `BROADNESS_PhaseSpaceStatistics` | Missing | Phase-space tests, metrics, and corrections. |
| `BROADNESS_PhaseSpace_RQA` | Missing | Recurrence metrics and plots. |
| `BROADNESS_Plot_ActivationTimeseries` | Missing | Statistical time-series plots. |
| `BROADNESS_RepeatedMeasuresANOVA` | Missing | Repeated-measures condition test. |
| `BROADNESS_SignificantSamplesToWindows` | Missing | Convert masks into time windows. |
| `BROADNESS_SpatialActivationClustering` | Missing | Clustering, results, plots, and optional NIfTI export. |
| `BROADNESS_SpatialPatternGroupStatistics` | Missing | Independent-group spatial inference. |
| `BROADNESS_SpatialPatternStatistics` | Missing | Within-participant spatial inference. |
| `BROADNESS_TemporalClusterCorrection` | Missing | Paired temporal cluster permutations. |
| `BROADNESS_TemporalClusterCorrection_ANOVA` | Missing | Repeated-measures temporal cluster permutations. |
| `BROADNESS_TemporalClusterCorrection_GroupANOVA` | Missing | Group Welch ANOVA with wild bootstrap. |
| `BROADNESS_TemporalClusterCorrection_Groups` | Missing | Pairwise group Welch tests with wild bootstrap. |
| `BROADNESS_Visualizer` | Existing, needs update | New figure mode, layout, export options, and output structure. |

The three requested examples are all pending: `BROADNESS_Example`,
`BROADNESS_Statistics_WithinParticipants_Example`, and
`BROADNESS_Statistics_MultipleGroupsParticipants_Example`.

## External functionality

| MATLAB external component | Python approach |
| --- | --- |
| `jadeR.m` | Implement or select an ICA method, document algorithm differences, and compare results on small fixtures. |
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
