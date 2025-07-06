function BROADNESS_Startup(path_home)

% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  STARTUP
% ========================================================================
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2024).
%  BROADband brain Network Estimation via Source Separation (BROAD-NESS). bioRxiv, 2024-10.
%  https://doi.org/10.1101/2024.10.31.621257
%
% ========================================================================
%
%  This function initializes the BROADNESS toolbox.
%
%  The user must provide `path_home` as a string, keeping the original 
%  subfolder structure of the BROADNESS toolbox unchanged.
%
%  Example usage:
%       path_home = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/CodeData/MIBSummerSchool2025';
%       addpath(path_home)
%       BROADNESS_Startup(path_home);
%
%  The function will add the relevant subfolders to the user's current MATLAB path.


% ------------------------------------------------------------------------
%  AUTHORS:
%  Leonardo Bonetti & Mattia Rosso
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  mattia.rosso@clin.au.dk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 05/07/2025

% ========================================================================



% Add main toolbox directories and NIFTI tools subfolder to MATLAB path
addpath(path_home);
addpath([path_home '/BROADNESS_Functions/'])
addpath([path_home '/BROADNESS_External/'])
addpath([path_home '/BROADNESS_External/NIfTI_20140122']);

end
