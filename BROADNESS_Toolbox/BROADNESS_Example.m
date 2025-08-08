%
% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2024).
%  BROADband brain Network Estimation via Source Separation (BROAD-NESS). bioRxiv, 2024-10.
%  https://doi.org/10.1101/2024.10.31.621257
%
%
% ========================================================================
%
%  This script demonstrates the application of BROADNESS 
%  (Broadband Brain Network Estimation via Source Separation).
%
%  The example dataset consists of data averaged across participants recorded
%  during an auditory (musical) memory recognition task.
%  The dataset in this format is currently not available online, while its previous
%  version (consisting of preprocessed epoched data) is avaiable here:
%  https://zenodo.org/records/10715160
%  The dataset in this format will be made publicly available soon.
%
%
%  PLEASE, NOTE:
%  The BROADNESS toolbox (PCA or ICA) can be run on data averaged across
%  participants to estimate the *spatial patterns* (brain networks). 
%  However, statistical analyses should NOT be performed on the group-averaged 
%  time series — instead, they should be computed at the participant level.
%
%  Recommended workflow:
%    1) Run BROADNESS_NetworkEstimation (PCA) or 
%       BROADNESS_AlternativeNetworkEstimation_ICA (ICA) 
%       on data averaged across participants to estimate the networks.
%
%    2) From the output, extract:
%          ActivationPatterns_BrainNetworks
%       (matrix of spatial patterns: [voxels × networks])
%
%    3) For each participant, project their own data onto these spatial patterns:
%          participant_TS = participant_data(:, 1:size(ActivationPatterns_BrainNetworks,2))' * ActivationPatterns_BrainNetworks;
%       This yields participant-specific network time series using the group-defined networks.
%
%    4) Perform your preferred statistical analyses on these participant-level
%       time series according to your experimental design.
%
%  This approach provides a fast, reliable way to estimate networks from 
%  group data, while ensuring that statistics are computed at the individual level.
%
%  We will soon release an example pipeline demonstrating these steps in code.
%  Until then, feel free to reach out for guidance:
%    Leonardo Bonetti  – leonardo.bonetti@clin.au.dk | leonardo.bonetti@psych.ox.ac.uk
%    Mattia Rosso      – mattia.rosso@clin.au.dk
%
%
%  The experimental task is described in detail both in the BROADNESS paper
%  referenced above and in the following paper:
%  Bonetti, L., Fernández-Rubio, G., Carlomagno, F., Dietz, M., Pantazis, D., Vuust, P.,
%  & Kringelbach, M. L. (2024). Spatiotemporal brain hierarchies of auditory memory
%  recognition and predictive coding. Nature Communications, 15(1), 4313.
%
%  The script loads the example datasets into the MATLAB workspace and 
%  calls three core functions:
%
% ------------------------------------------------------------------------
%  FUNCTIONS OVERVIEW:
% ------------------------------------------------------------------------
%  - 0) BROADNESS_Startup() 
%       Initializes the environment.
%
%  - 1) BROADNESS_NetworkEstimation() 
%       Performs PCA on the event-related field/potentials to derive the
%       underlying brain networks.
%
%  - 2) BROADNESS_Visualizer() 
%       Takes as input selected outputs from `BROADNESS_NetworkEstimation` 
%       to visualize a number of features such as brain network time series
%       and topographies.
%
%  - 3) BROADNESS_PhaseSpace_RQA()
%       Computes phase space embedding and RQA on brain networks time
%       series [it works for both the time series extracted by
%       BROADNESS_NetworkEstimation (PCA) and BROADNESS_AlternativeNetworkEstimation_ICA (ICA)]
%
%  - 4) BROADNESS_SpatialGradients()
%       Computes spatial gradients embedding and clustering on the spatial
%       activation patterns of the brain networks
%
%  - 5) BROADNESS_AlternativeNetworkEstimation_ICA()
%       Alternative computation of brain networks using ICA
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Leonardo Bonetti, Chiara Malvaso & Mattia Rosso
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  chiara.malvaso@studio.unibo.it
%  mattia.rosso@clin.au.dk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Department of Physics, University of Bologna
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 08/08/2025
% ========================================================================



%% 0) STARTUP

% Simply download the BROADNESS Toolbox folder and place it in your working directory,
% making sure not to alter the structure of its functions, subfolders, or files.

clear 
close all
clc

% Setup directories
path_home = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/CodeData/MIBSummerSchool2025/BROADNESS_Toolbox';
addpath(path_home)
BROADNESS_Startup(path_home);

%%

%% 1) PERFORM BROADNESS (ONLY ESSENTIAL INPUTS)

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply loading data (brain_voxel-by-time matrix) and time
load([path_home '/Data_Example.mat'])  % Resting-state data

%%% ------------------ COMPUTATION --------------------- %%%

% Run BROADNESS network estimation (default parameters)
BROADNESS = BROADNESS_NetworkEstimation(data, time);

%%

%% 1b) PERFORM BROADNESS (ALTERNATIVE SCENARIO WITH OPTIONAL INPUTS)

% This section demonstrates the same function as above,  
% but with optional settings provided. Any missing arguments  
% will automatically use their default values.  
% NOTE: YOU ONLY NEED TO RUN ONE SECTION: EITHER THIS ONE OR THE PREVIOUS ONE.  


%%% ------------------- USER SETTINGS ------------------- %%%

% Simply loading data (brain_voxel-by-time matrix) and time
load([path_home '/Data_Example.mat'])  % Resting-state data
% Optional arguments
time_window = [0.350 1.750];
permutations_num = 10;
randomization = 1;
sign_eigenvect = '0';

%%% ------------------ COMPUTATION --------------------- %%%

BROADNESS = BROADNESS_NetworkEstimation(data, time, ...
                                'time_window', time_window, 'permutations_num', permutations_num, 'randomization', randomization, 'sign_eigenvect', sign_eigenvect); %% Call with optional parameters

%%

%% 2) BROADNESS VISUALIZATION

%This section generates 5 plots, described as follows:
%   #1) Dynamic brain activity map of the original data
%   #2) Variance explained by the networks
%   #3) Time series of the networks
%   #4) Activation patterns of the networks (3D)
%   #5) Activation patterns of the networks (nifti images)


%%% ------------------- USER SETTINGS ------------------- %%%

% Minimal user settings: output folder and MNI coordinates of original brain voxel data
Options = [];
Options.name_nii = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/CodeData/MIBSummerSchool2025'; %output folder
load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;

%%% ------------------ COMPUTATION --------------------- %%%

% Visualize brain networks features
BROADNESS_Visualizer(BROADNESS,Options)

%%

%% 2b) BROADNESS VISUALIZATION (ALTERNATIVE SCENARIO WITH OPTIONAL INPUTS)

% This section demonstrates the same function as above,  
% but with optional settings provided. Any missing arguments  
% will automatically use their default values.  
% NOTE: YOU ONLY NEED TO RUN ONE SECTION: EITHER THIS ONE OR THE PREVIOUS ONE.  

%%% ------------------- USER SETTINGS ------------------- %%%

% Minimal user settings: output folder
Options = [];
Options.name_nii = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/CodeData/MIBSummerSchool2025'; %output folder
load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;
Options.WhichPlots = [0 0 1 0 0]; %which plots to be generated
Options.ncomps = [1:3]; %indices of PCs to be plotted (all plots)
Options.ncomps_var = 60; %number of PCs to be plotted (only in Variance plot)
Options.Labels = {'Memorized','NewT1','NewT2','NewT3','NewT4'}; %experimental condition labels
Options.color_PCs = [
    0.4,    0.761,  0.647;   % teal-green
    0.988,  0.553,  0.384;   % coral
    0.553,  0.627,  0.796;   % periwinkle
    0.906,  0.541,  0.765;   % pink-purple
    0.651,  0.847,  0.329;   % green
    1.000,  0.851,  0.184;   % yellow
    0.898,  0.769,  0.580;   % beige
    0.702,  0.702,  0.702    % gray
]; %RGB code colors for PCs
Options.color_conds = [
    0.106,  0.620,  0.467;   % green
    0.851,  0.373,  0.008;   % orange
    0.459,  0.439,  0.702;   % purple
    0.906,  0.161,  0.541;   % magenta
    0.400,  0.651,  0.118;   % lime green
    0.902,  0.671,  0.008;   % gold
    0.651,  0.463,  0.114;   % brown
    0.4,    0.4,    0.4      % gray
]; %RGB code colors for experimental conditions

% If you wish to remove the cerebellum voxels (not included in 3D brain template (#4)), please set 'remove_cerebellum_label' to 1
% NOTE: This removal works only for 8mm brain
remove_cerebellum_label = 0;
if remove_cerebellum_label == 1
    load([path_home '/BROADNESS_External/cerebellum_coords.mat']); %only cerebellar voxels
    % Remove cerebellar voxels since they are not included in the 3D brain template (#4)
    [~, idx_cerebellum] = ismember(MNI8, cerebellum_coords, 'rows');  % find cerebellum indexes in MNI coordinates matrix (all voxels)
    MNI8(idx_cerebellum~=0,:) = nan; %assigning nans to MNI coordinates matrix
    Options.MNI_coords = MNI8; %assigning the MNI coordinates of your data for visualization purposes (both 3D main template (#4) and nifti images (#5))
end

%%% ------------------ COMPUTATION --------------------- %%%

% Visualize brain networks features
BROADNESS_Visualizer(BROADNESS,Options)

%%

%% *** ANALYSIS ON SPATIAL-TEMPORAL FEATURES OF BROADNESS BRAIN NETWORKS ***

%%

%% 3) PHASE SPACE EMBEDDING AND RECURRENCE QUANTIFICATION ANALYSIS 

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply use the structure outputted by the BROADNESS_NetworkEstimation function
% Additional optional inputs can be provided, as described in the function. 

RQA_BROADNESS = BROADNESS_PhaseSpace_RQA(BROADNESS,'principalcomps',[1:2],'threshold',0.1,'video','on','figure','on');

%%

%% 4) SPATIAL GRADIENTS EMBEDDING AND CLUSTERING ANALYSIS

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply use the structure outputted by the BROADNESS_NetworkEstimation function
% Additional optional inputs can be provided, as described in the function. 

SPATIAL_GRADIENTS_BROADNESS = BROADNESS_SpatialGradients(BROADNESS,'principalcomps',[1:2],'evalclusters',1,'scatterplots','all','outpath','/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/CodeData/MIBSummerSchool2025');

%%

%% *** BROADNESS ALTERNATIVE NETWORK ESTIMATION - INDEPENDENT COMPONENT ANALYSIS (ICA) ***

%% 5) BROADNESS ICA COMPUTATION

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply loading data (brain_voxel-by-time matrix) and time
load([path_home '/Data_Example.mat'])  % Resting-state data
% Additional optional inputs can be provided, as described in the function. 

BROADNESS_ICA = BROADNESS_AlternativeNetworkEstimation_ICA(data, time, 'icacomps', 5, 'total_varexp', 95);

%% 2) BROADNESS VISUALIZATION (USED FOR ICA)

%This section generates 5 plots, described as follows:
%   #1) Dynamic brain activity map of the original data
%   #2) Variance explained by the networks
%   #3) Time series of the networks
%   #4) Activation patterns of the networks (3D)
%   #5) Activation patterns of the networks (nifti images)


%%% ------------------- USER SETTINGS ------------------- %%%

% Minimal user settings: output folder
Options = [];
Options.name_nii = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/CodeData/MIBSummerSchool2025'; %output folder
load([path_home '/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;
Options.WhichPlots = [0 0 1 0 0]; %which plots to be generated
Options.ncomps = [3]; %indices of PCs to be plotted (all plots)
Options.Labels = {'Memorized','NewT1','NewT2','NewT3','NewT4'}; %experimental condition labels


%%% ------------------ COMPUTATION --------------------- %%%

% Visualize brain networks features
BROADNESS_Visualizer(BROADNESS_ICA,Options)

%%

%%

%%

%% FINAL REMARKS

% This demonstration uses data from averaged participants.  
% If you want to compute statistics, the best solution is to use the
% activation patterns computed here (in PCA on the average across participants)
% to calculate the time series independently for each participant (i.e.
% multiplying the "grup-level" activation patterns by the data,
% independnely for each participant).

% Please check the BROADNESS GitHub repository for new releases.  
% https://github.com/leonardob92/BROADNESS_MEG_AuditoryRecognition.git
% Feel free to reach out to us if you need guidance or consultation.  
% Leonardo Bonetti: leonardo.bonetti@clin.au.dk
%                   leonardo.bonetti@psych.ox.ac.uk
% Mattia Rosso:     mattia.rosso@clin.au.dk
% Chiara Malvaso:   chiara.malvaso@studio.unibo.it
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2024).
%  BROADband brain Network Estimation via Source Separation (BROAD-NESS). bioRxiv, 2024-10.
%  https://doi.org/10.1101/2024.10.31.621257

%%
