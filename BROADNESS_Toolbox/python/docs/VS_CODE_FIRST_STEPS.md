# BROADNESS in VS Code: first steps

This guide starts with a fresh computer and ends with a small working BROADNESS
analysis. It uses a `.py` script with sections, like a MATLAB live script. The
example uses generated data, so you do **not** need the large research dataset.

**Already ran setup successfully on this computer?** Skip section 2. Open the
`python` folder in VS Code, select its `.venv` interpreter in section 3, and
run the four cells in section 4. You do not need to install Python again.

## 1. Get the software and repository

1. Install [Visual Studio Code](https://code.visualstudio.com/download).
2. Download the BROADNESS repository using **Code → Download ZIP** on its
   [GitHub page](https://github.com/leonardob92/BROADNESS_MEG_AuditoryRecognition),
   then extract it. If you already cloned the repository with GitHub Desktop,
   use that existing folder instead.
3. In VS Code, choose **File → Open Folder** and open
   `BROADNESS_MEG_AuditoryRecognition/BROADNESS_Toolbox/python`.
4. Check the two required extensions yourself:
   - Press **Cmd+Shift+X** on a Mac or **Ctrl+Shift+X** on Windows/Linux. This
     opens **Extensions** (the four-squares icon on the left).
   - Search for **`ms-python.python`**. Open the result named **Python** by
     **Microsoft**. If you see **Install**, click it. If you see **Enable** or
     **Enable (Workspace)**, click that. If you see **Disable** or **Uninstall**,
     it is already installed and enabled.
   - Clear the search box and search for **`ms-toolsai.jupyter`**. Open
     **Jupyter** by **Microsoft** and make the same check.
   - If VS Code asks you to **Reload Window** or **Restart Extensions**, do so.
     Then reopen this folder if needed. Both extension pages should now offer
     **Disable** or **Uninstall**, rather than **Install** or **Enable**.

This repository recommends those extensions, but a recommendation does not
install them automatically. The check above confirms their actual state.

VS Code is the editor. Python is the program that runs the analysis, and
`.venv` will be a private Python environment for this toolbox.

## 2. Check Python and install BROADNESS

In VS Code, choose **Terminal → New Terminal**. The prompt should be in the
`python` folder you opened above. Run the command for your operating system:

| Computer | Command to type in the VS Code terminal |
| --- | --- |
| macOS or Linux | `bash setup_vscode.sh` |
| Windows (PowerShell or Command Prompt) | `.\setup_vscode.bat` |

The setup script checks for **Python 3.11 through 3.14**, creates `.venv`, and
installs BROADNESS, plotting support, and the VS Code kernel. A successful run
ends with **`BROADNESS and the VS Code kernel are ready`**. It does not change
your other Python projects.

If the script says **No compatible Python was found**:

1. Open the [Python 3.14 installer page](https://www.python.org/downloads/release/python-31410/)
   and download the installer for your operating system. Do not select Python
   3.15 for this toolbox yet.
2. On a Mac, double-click the downloaded `.pkg` and follow the installer
   screens. On Windows, run the installer and enable **Add Python to PATH** if
   it offers that choice.
3. Close the VS Code terminal, open **Terminal → New Terminal**, and type the
   same setup command from the table above. Wait for the success message.

System Python installation may ask for your operating-system password. The
shell or batch script detects Python before trying to run Python code.

If a previous installation was interrupted, run the same setup command again.
An existing compatible `.venv` is reused. An incompatible one is moved to a
`.venv-previous-*` backup before setup creates a replacement. If setup reports
a network error, connect to the internet and retry. Installation downloads
Python packages once; the quickstart uses very little memory and no participant
data.

**If you installed Python 3.15:** the first setup script could try to compile
`h5py` and fail because a compatible prebuilt wheel was unavailable. You do
not need to uninstall Python 3.15. Install Python 3.14 alongside it, or point
the macOS/Linux setup script at an existing compatible Python with
`BROADNESS_PYTHON=/full/path/to/python bash setup_vscode.sh`. The script will
move the incomplete 3.15 `.venv` aside and build a new one.

## 3. Tell VS Code which Python to use

1. Press **Cmd+Shift+P** on macOS or **Ctrl+Shift+P** on Windows/Linux.
2. Type **Python: Select Interpreter** and press Enter.
3. Choose the interpreter inside this folder's `.venv`. On macOS/Linux its
   path ends in `.venv/bin/python`; on Windows it ends in
   `.venv\Scripts\python.exe`.
4. Look at the lower right of VS Code: the selected Python interpreter should
   say `.venv`. If it does not appear, use **Enter interpreter path** and select
   that file directly.

On a Mac, the selected path should be inside the folder you opened in step 1,
ending in `BROADNESS_Toolbox/python/.venv/bin/python`. The global `python3`
command can point to a different Python; the selected `.venv` path is what
matters for BROADNESS.

You do not need to activate `.venv` manually to use **Run Cell**. The VS Code
Jupyter extension uses the selected environment as the kernel. If VS Code asks
you to **Select Kernel**, choose **Python Environments → .venv**.

## 4. Run the sections, one at a time

1. Open [`examples/00_vscode_quickstart.py`](../examples/00_vscode_quickstart.py)
   in the VS Code file explorer.
2. Above `# %% 1. Import tools`, click **Run Cell**. The **Python Interactive**
   window opens. Wait for the `BROADNESS is ready` message.
3. Click **Run Cell** above section 2, then section 3, then section 4. Keep the
   same Interactive Window open. The variables produced by earlier sections
   remain available to later sections.
4. Section 4 makes a simple plot and prints where to inspect the variables.

You have finished the first lesson when you see a data shape of `(12, 80, 2, 6)`,
a network time-series shape of `(80, 4, 2, 6)`, and a two-line plot. If a cell
reports that a name such as `data` or `pca` is missing, restart the Interactive
Window and run cells 1–4 in order again.

If **Run Cell** is missing, check that the file name ends in `.py` and that the
Python and Jupyter extensions are installed. You can also open the Command
Palette and run **Jupyter: Run Current File in Python Interactive Window** to
execute the whole quickstart.

## 5. Find the MATLAB-style workspace

In the **Python Interactive** window, click **Variables** in its toolbar. This
opens the **Variables explorer**. You should see `data`, `time`, `pca`,
`component_1`, and other names. `data` has shape `(12, 80, 2, 6)`: 12 sources,
80 time samples, 2 conditions, and 6 participants. Double-click a NumPy array
to inspect its values in the **Data Viewer**. `pca` is a result object; use
`pca.TimeSeries_BrainNetworks` or `pca.Variance_BrainNetworks` in a later cell
to inspect its fields.

The Variables explorer shows the **current kernel's memory**, like MATLAB's
workspace for a running session. Closing or restarting that kernel clears the
variables. After changing an earlier section, run that section and all sections
that depend on it again. For a clean start, choose **Restart** in the
Interactive Window toolbar and run cells 1–4 in order.

## 6. Download real example data (optional)

The first lesson above needs no download. To try the published data, open the
[BROAD-NESS dataset on Zenodo](https://zenodo.org/records/17048137). Create a
folder called `Data` at the **top level of the repository**, alongside
`BROADNESS_Toolbox`, and put the extracted `.mat` files directly in it:

```text
BROADNESS_MEG_AuditoryRecognition/
├── BROADNESS_Toolbox/
│   ├── matlab/
│   └── python/
└── Data/
    ├── DataReduced_AveragedOverParticipants_Example.mat
    ├── SUBJ_0002.mat
    ├── SUBJ_0007.mat
    └── ...
```

Choose the files for the kind of analysis you want:

| Goal | Download from Zenodo | Place in `Data/` |
| --- | --- | --- |
| Try group-average PCA and plots | `DataReduced_AveragedOverParticipants_Example.mat` (about 51 MB) | That one `.mat` file |
| Try participant-level loading and PCA | `Dataset_1_IndividualParticipants.zip` or `Dataset_2_IndividualParticipants.zip`, plus the averaged `.mat` file for its `time` vector | Extract a few `SUBJ_*.mat` files and the averaged `.mat` file directly into `Data/` |

Each participant file is about 140 MB. Zenodo supplies them **inside large ZIP
archives** (about 6.0 GB and 5.6 GB), so even if you extract only a few
participants, you may need to download an entire archive first. Do not put
the ZIP itself, or an extra nested folder, where the examples expect the
`.mat` files. Participant-level **statistics** need enough participants for
the chosen comparison; a few files are useful for learning the loading and
PCA workflow, not for reproducing the paper's inference.

The three Python examples follow the corresponding MATLAB examples. Start with
the [general tutorial](tutorials/01_GENERAL.md), then the
[within-participant tutorial](tutorials/02_WITHIN_PARTICIPANTS.md) or
[multiple-group tutorial](tutorials/03_MULTIPLE_GROUPS.md). Their default data
location is the `Data/` folder shown above. If you store the files elsewhere,
change the `data_path` setting in section 0 of each script to your folder, for example
`data_path = Path("/absolute/path/to/my/data")` after `from pathlib import Path`.
On Windows, use a raw string such as `Path(r"C:\my_data\BROADNESS")`.

For BROADNESS, the main input order is **sources × time × conditions ×
participants**. Public component numbers start at 1, like MATLAB; when you
directly index a NumPy array, position 0 is the first item. Keep these data
files on your computer: the repository ignores `Data/`, so they will not
appear in GitHub Desktop commits.

You can run these `.py` sections in VS Code without creating a Jupyter notebook.
Each script remains an ordinary Python file that can also run from top to bottom.

VS Code references: [Python setup](https://code.visualstudio.com/docs/python/python-tutorial),
[Interactive Window and Variables explorer](https://code.visualstudio.com/docs/python/jupyter-support-py),
and [kernel selection](https://code.visualstudio.com/docs/datascience/jupyter-kernel-management).
