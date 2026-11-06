# BROADNESS Toolbox

Choose the implementation for your environment:

- [`matlab/`](matlab/) contains `BROADNESS_Startup.m`, all MATLAB functions,
  required external files, and the general and two statistics examples.
- [`python/`](python/) contains the installable Python package, the
  [VS Code first-steps guide](python/docs/VS_CODE_FIRST_STEPS.md), and three
  [real-data tutorials](python/README.md). See the
  [porting checklist](python/docs/PORTING_STATUS.md) for verification status.

The examples look for research data in the repository's `Data/` directory; the
[Python first-steps guide](python/docs/VS_CODE_FIRST_STEPS.md#6-download-real-example-data-optional)
explains how to get the [published dataset](https://zenodo.org/records/17048137)
and place the files there.
Generated figures and analysis files go into `BROADNESS_Toolbox/Output/`, which
is local and ignored by Git.
