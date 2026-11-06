# BROADNESS Toolbox

Choose the implementation for your environment:

- [`matlab/`](matlab/) contains `BROADNESS_Startup.m`, all MATLAB functions,
  required external files, and the general and two statistics examples.
- [`python/`](python/) contains the installable Python package, tests, and
  [porting checklist](python/docs/PORTING_STATUS.md). The Python port is still
  in progress.

The examples look for research data in the repository's `Data/` directory.
Generated figures and analysis files go into `BROADNESS_Toolbox/Output/`, which
is local and ignored by Git.
