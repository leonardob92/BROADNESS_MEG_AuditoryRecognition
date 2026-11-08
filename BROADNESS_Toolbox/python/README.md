# BROADNESS for Python

The Python toolbox has the same 20 analysis functions as the MATLAB toolbox.
The examples use `# %%` sections, so you can run them one section at a time in
Visual Studio Code and inspect their variables in the Interactive Window.

## Start here

1. Follow [VS Code first steps](docs/VS_CODE_FIRST_STEPS.md) to install the
   toolbox, select its Python interpreter, and run a tiny synthetic-data lesson.
2. Download the [published data from Zenodo](https://zenodo.org/records/17048137)
   and place the `.mat` files in the repository's top-level `Data/` folder as
   shown in the [data-folder guide](docs/VS_CODE_FIRST_STEPS.md#6-download-real-example-data-optional).
3. Choose a real-data tutorial:

| Tutorial | Script | Data needed |
| --- | --- | --- |
| [1. General networks](docs/tutorials/01_GENERAL.md) | [`01_BROADNESS_Example.py`](examples/01_BROADNESS_Example.py) | Averaged `.mat` file; participants optional |
| [2. Conditions within participants](docs/tutorials/02_WITHIN_PARTICIPANTS.md) | [`02_BROADNESS_Statistics_WithinParticipants_Example.py`](examples/02_BROADNESS_Statistics_WithinParticipants_Example.py) | Averaged file plus at least four `SUBJ_*.mat` files |
| [3. Independent groups](docs/tutorials/03_MULTIPLE_GROUPS.md) | [`03_BROADNESS_Statistics_MultipleGroupsParticipants_Example.py`](examples/03_BROADNESS_Statistics_MultipleGroupsParticipants_Example.py) | Averaged file plus at least six `SUBJ_*.mat` files |

The general example uses the full averaged dataset. The two statistics scripts
start in **tutorial mode**, sampling sources and time points to keep the first
run manageable. Those subsets and the arbitrary groups in tutorial 3 are for
learning the software, not scientific inference.

## Further information

- [Function and installation guide](docs/FUNCTION_GUIDE.md): inputs, outputs,
  plotting, optional dependencies, indexing, and scientific conventions.
- [Porting status](docs/PORTING_STATUS.md): MATLAB correspondence and remaining
  numerical comparison work.
- [MATLAB toolbox](../matlab/): original functions and example scripts.

Participant files and generated outputs remain local and are ignored by Git.
