function [ED] = BROADNESS_EffectiveDimensionality(eigenspectrum)
%%
% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  EFFECTIVE DIMENSIONALITY ESTIMATION
% ========================================================================
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2025). Advanced Science. 
%  BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in Auditory Memory Networks.
%  https://doi.org/10.1002/advs.202507878
%
% ========================================================================
%
%  This function returns the effective dimensionality of the decomposition,
%  providing the number of principal components (brain networks) that 
%  contain relevant information and that is recommended to use for additional
%  analysis and to report in the results section of your study.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - eigenspectrum:     Vector with variance explained by the different PC (brain networks).
%                       It corresponds to: BROADNESS.Variance_BrainNetworks
%                       Obtained from the BROADNESS_NetworkEstimation function.
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - ED:                Effective dimensionality (i.e., number of PC to be kept)
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso, Chiara Malvaso & Leonardo Bonetti 
%  mattia.rosso@clin.au.dk
%  chiara.malvaso@studio.unibo.it
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Department of Physics, University of Bologna
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 23/08/2025
%
% ========================================================================
%
%%





%% Compute effective dimensionality (ED)

num   = sum(eigenspectrum, 1).^2;       % [1 x nFrex x nSubs]
denom = sum(eigenspectrum.^2, 1);       % [1 x nFrex x nSubs]

% Guard against division by zero
denom = max(denom, realmin);

% Compute (sum(lambda))^2 / sum(lambda.^2) along component dimension
ED = round(squeeze(num ./ denom));             % [nFrex x nSubs]


end

