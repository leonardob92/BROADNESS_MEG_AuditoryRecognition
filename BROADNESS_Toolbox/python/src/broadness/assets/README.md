# BROADNESS anatomical visualization asset

`MNI152_T1_1mm_brain_surface.npz` contains a reduced triangular surface in MNI
world coordinates. It was generated once from the project-provided
`MNI152_T1_1mm_brain.nii.gz` brain mask and is used only as the neutral
anatomical background of the Python spatial-pattern visualization.

The archive contains 6,013 vertices and 12,000 triangular faces. It contains no
participant data, network estimates, or experimental measurements. Network
activation patterns remain at their original source coordinates and resolution.

The small 8-mm MNI152 NIfTI templates and MATLAB coordinate files are also
bundled here so Python visualization and future spatial analyses can work after
installation, without locating the MATLAB toolbox. These are anatomical
reference assets, not participant data.
