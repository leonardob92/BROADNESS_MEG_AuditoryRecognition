BROADband brain Network Estimation via Source Separation (BROAD-NESS)

Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno, F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2025). Advanced Science, e07878. 
BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in Auditory Memory Networks.
https://doi.org/10.1002/advs.202507878

Additional relevant codes and functions are available here: https://github.com/leonardob92/LBPD-1.0.git

Abstract: Understanding perceptual and cognitive processes typically relies on event-related designs, where tasks are repeated in brief time windows. While it is established that the brain is organised into functional networks, these are often assessed relying on predefined regions of interest and sophisticated techniques originally designed for resting-state data. This hinders interpretation of event-related networks and does not fully capitalise on the high temporal and spatial resolution of magnetoencephalography (MEG). Here, we introduce BROADband brain Network Estimation via Source Separation (BROAD-NESS), a novel method that uses principal component analysis (PCA) on source-reconstructed MEG data at the voxel level to identify meaningful brain networks. We applied BROAD-NESS to MEG data from 83 participants engaged in a long-term musical sequence recognition task. We identified two main networks accounting for 88% of the variance (72% and 16%, respectively), both including the auditory cortices. The first network also involved the medial cingulate gyrus and was linked to early and sustained auditory processes. The second network, in addition to the auditory cortices, comprised the prefrontal and hippocampal regions, inferior temporal cortex and insula. This network’s time series showed a slow positive peak after each tone, suggesting a process of matching sounds to memory traces of previously heard melodies, alongside a faster negative peak indicating prediction errors when the match was incorrect. In conclusion, BROAD-NESS effectively identified simultaneous, broadband brain networks involved in long-term musical sequence recognition, enhancing our understanding of memory-related brain dynamics and offering a novel tool to investigate brain networks in event-related designs.


Please, note #1:

-BroadNess_APR2020.m contains code for the preprocessing, source reconstruction and PCA computed on the data averaged across participants.

-BroadNess_APR2020_2.m contains statistical testing and several tests performed in relation to the PCA computation.

-If you wish to use our code for running PCA and derive brain networks, please refer to the function: PCA_LBPD.mat which is available here: https://github.com/leonardob92/LBPD-1.0 and to the following short script showing how to use it: "PCA_UseMe.m"

Please, note #2:

-The BROADNESS Toolbox has separate [MATLAB](BROADNESS_Toolbox/matlab) and [Python](BROADNESS_Toolbox/python) folders. The MATLAB functions and examples are in `BROADNESS_Toolbox/matlab`; the Python package is being updated to cover the same functionality.
-The BROADNESS Toolbox has separate [MATLAB](BROADNESS_Toolbox/matlab) and [Python](BROADNESS_Toolbox/python) folders. The MATLAB functions and examples are in `BROADNESS_Toolbox/matlab`. New Python users can start with the [VS Code first-steps guide](BROADNESS_Toolbox/python/docs/VS_CODE_FIRST_STEPS.md) and then choose one of the [three real-data tutorials](BROADNESS_Toolbox/python/README.md).

The conceptualization and code implementation of the Toolbox were carried out in collaboration with Dr. Mattia Rosso (https://github.com/mattiaRosso92).
Corresponding authors: Leonardo Bonetti - leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk; Mattia Rosso - mattia.rosso@clin.au.dk
