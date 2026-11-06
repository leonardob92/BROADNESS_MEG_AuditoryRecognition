# BROADNESS for MATLAB

This folder contains the MATLAB toolbox, its required external files, and the
three examples. To use its functions from another MATLAB script, set the path
to this folder and call startup:

```matlab
toolbox_path = '/path/to/BROADNESS_Toolbox/matlab';
addpath(toolbox_path)
BROADNESS_Startup(toolbox_path)
```

Start with `BROADNESS_Example.m` for the general workflow, or use
`BROADNESS_Statistics_WithinParticipants_Example.m` and
`BROADNESS_Statistics_MultipleGroupsParticipants_Example.m` for the two
statistics workflows. The examples expect input files in the repository's
`Data/` folder; adjust their `data_path` if yours are elsewhere. Generated
results are saved under `BROADNESS_Toolbox/Output/` and are ignored by Git.
