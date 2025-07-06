%
% ========================================================================
%  BROADBAND (BRAIN) NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
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
%  - BROADNESS_Startup(path_home) 
%    Initializes the environment. Given `path_home` as input, it sets up 
%    the necessary directories.
%
%  - BROADNESS_NetworkEstimation(...) 
%    Performs PCA on the event-related field/potentials to derive the
%    underlying brain networks.
%
%  - BROADNESS_Visualizer(...) 
%    Takes as input selected outputs from `BROADNESS_NetworkEstimation` 
%    to visualize a number of features such as brain network time series
%    and topographies.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Leonardo Bonetti & Mattia Rosso
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  mattia.rosso@clin.au.dk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 05/07/2025
% ========================================================================



%% STARTUP

clear 
close all
clc

% Setup directories
path_home = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/CodeData/MIBSummerSchool2025';
addpath(path_home)
BROADNESS_Startup(path_home);

%%

%% PERFORM BROADNESS (ONLY ESSENTIAL INPUTS)

%%% ------------------- USER SETTINGS ------------------- %%%

% Simply loading data (brain_voxel-by-time matrix) and time
load([path_home '/Data_Example.mat'])  % Resting-state data

%%% ------------------ COMPUTATION --------------------- %%%

% Run BROADNESS network estimation (default parameters)
BROADNESS = BROADNESS_NetworkEstimation(data, time);

%%

%% PERFORM BROADNESS (ALTERNATIVE SCENARIO WITH OPTIONAL INPUTS)

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

%% BROADNESS VISUALIZATION

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
load([path_home '/BROADNESS_Toolbox/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;

%%% ------------------ COMPUTATION --------------------- %%%

% Visualize brain networks features
BROADNESS_Visualizer(BROADNESS,Options)

%%

%% BROADNESS VISUALIZATION (ALTERNATIVE SCENARIO WITH OPTIONAL INPUTS)

% This section demonstrates the same function as above,  
% but with optional settings provided. Any missing arguments  
% will automatically use their default values.  
% NOTE: YOU ONLY NEED TO RUN ONE SECTION: EITHER THIS ONE OR THE PREVIOUS ONE.  

%%% ------------------- USER SETTINGS ------------------- %%%

% Minimal user settings: output folder
Options = [];
Options.name_nii = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/CodeData/MIBSummerSchool2025'; %output folder
load([path_home '/BROADNESS_Toolbox/BROADNESS_External/MNI152_8mm_coord_dyi.mat']); %all voxels MNI coordinates
Options.MNI_coords = MNI8;
Options.WhichPlots = [0 0 0 1 0]; %which plots to be generated
Options.ncomps = 5; %number of PCs to be plotted (all plots)
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
    load([path_home '/BROADNESS_Toolbox/BROADNESS_External/cerebellum_coords.mat']); %only cerebellar voxels
    % Remove cerebellar voxels since they are not included in the 3D brain template (#4)
    [~, idx_cerebellum] = ismember(MNI8, cerebellum_coords, 'rows');  % find cerebellum indexes in MNI coordinates matrix (all voxels)
    MNI8(idx_cerebellum~=0,:) = nan; %assigning nans to MNI coordinates matrix
    Options.MNI_coords = MNI8; %assigning the MNI coordinates of your data for visualization purposes (both 3D main template (#4) and nifti images (#5))
end

%%% ------------------ COMPUTATION --------------------- %%%

% Visualize brain networks features
BROADNESS_Visualizer(BROADNESS,Options)

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
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2024).
%  BROADband brain Network Estimation via Source Separation (BROAD-NESS). bioRxiv, 2024-10.
%  https://doi.org/10.1101/2024.10.31.621257

%%
