function [BROADNESS_ICA] = BROADNESS_AlternativeNetworkEstimation_ICA(data, time, varargin)
%%
% =========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  ICA-BASED BRAIN NETWORK ESTIMATION
% =========================================================================
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2025). Advanced Science. 
%  BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in Auditory Memory Networks.
%  https://doi.org/10.1002/advs.202507878
%
% =========================================================================
%
%  This function applies Independent Component Analysis (ICA) to extract
%  brain networks from a multivariate M/EEG dataset based on the statistical
%  independence of their time series.
%
%  The approach is optimized for source-reconstructed ERP/F data
%  Please, note, that this function can be used at the single trial level,
%  but we suggest to compute it using data averaged across trials.
%  Same concept applies to participants.
%
%  If the user does not specify the number of ICs, we default to the number
%  of PCA components explaining at least 95% of the variance (configurable).
%
% -------------------------------------------------------------------------
%  INPUTS:

%    data   :  DATA MATRIX, either:
%                   - 2D: sources × time-points  (or channels × time-points)
%                   - 3D: sources × time-points × conditions
%                   - 4D: sources × time-points × conditions × participants
%                   Note on PCA computation:
%                   - If 3D: average across conditions   → reduces to 2D
%                   - If 4D: average across participants and then across conditions → reduces to 2D
%                   General note:
%                   Depending on experimental needs, BROADNESS can be computed:
%                   - after averaging across conditions/participants, OR
%                   - independently for each condition and participant.
%                   This is a choice that the user should make, after thinking carefully.
%
%    time   : time vector (seconds), length must equal size(data,2)
%
%  Name-value pairs:
%    'icacomps'     : [] (default). If provided, use this integer #ICs.
%    'total_varexp' : 95% (default). Used to pick #ICs via PCA when
%                     'icacomps' is not given. Must be in [0,100].
%
%  Behavior:
%    If 'icacomps' is provided, 'total_varexp' is ignored.
%    Otherwise PCA is used to pick the smallest m with cumulative variance
%    >= total_varexp; ICA is then run with m components.
%
% -------------------------------------------------------------------------
%  OUTPUT (BROADNESS_ICA struct):
%    .TimeSeries_BrainNetworks      : IC time courses (IC × time × cond)
%    .ActivationPatterns_BrainNetworks : Spatial activation patterns
%                                        (voxels × IC). Computed as
%                                        (unmixing weights) * Cov(data), transposed.
%    .numICAcomps                   : number of ICs used (m)
%    .Time                          : time vector
%    .OriginalData                  : original input data
%
% -------------------------------------------------------------------------
%  AUTHORS:
%  Chiara Malvaso, Mattia Rosso & Leonardo Bonetti 
%  chiara.malvaso@studio.unibo.it
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Department of Physics, University of Bologna
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 08/08/2025
% =========================================================================
%
% NOTE: This function uses the JADE ICA algorithm implemented in the jadeR function.
% JADE was developed by Jean-François Cardoso.
% Cardoso, J.-F., & Souloumiac, A. (1993). Blind beamforming for non-Gaussian signals.
% IEE Proceedings F (Radar and Signal Processing), 140(6), 362–370. https://doi.org/10.1049/ip-f-2.1993.0054

%%




%% ----------------------------- Handle inputs -----------------------------
disp('Checking inputs')

opts = struct('total_varexp', 0.95, 'icacomps', []);
opts = parse_name_value_pairs(opts, varargin{:});
total_varexp = opts.total_varexp;
icacomps     = opts.icacomps;

if ~isnumeric(total_varexp) || ~isscalar(total_varexp) || total_varexp <= 0 || total_varexp > 100
    error('"total_varexp" must be a numeric scalar between 0 and 100');
end

total_varexp = total_varexp/100;

sz = size(data);
if ~(ndims(data) == 2 || ndims(data) == 3 || ndims(data) == 4)
    error('Data must be 2D (voxels × time) or 3D (voxels × time × conditions) or 4D (voxels × time × conditions x participants)');
end
if size(data,2) ~= numel(time)
    error('size(data,2) must match length(time)');
end

%% ------------------------- Preprocess data ------------------------------

% Average across conditions if 3D
if ndims(data) == 4
    data_temp = data; % store the data
    data = mean(data, 4); % average across participants
    data = mean(data, 3); % average across conditions
elseif ndims(data) == 3
    data_temp = data; % same data stored with a different name for avoiding issues later (barbaric yet effective solution)
    data = mean(data, 3); % average across conditions
end

% Demean across time for each voxel
data_demeaned = data - mean(data,2);

%% -------- Determine number of ICA components (m) ------------------------

if ~isempty(icacomps)
    % User-specified number of components
    if ~isscalar(icacomps) || ~isnumeric(icacomps) || icacomps <= 0 || fix(icacomps) ~= icacomps
        error('"icacomps" must be a positive integer scalar.');
    end
    maxAllowed = min(size(data_demeaned));
    if icacomps > maxAllowed
        error('"icacomps" (%d) exceeds allowable components (%d).', icacomps, maxAllowed);
    end
    disp('Computing ICA with the requested number of components')
    m = icacomps;
else
    % PCA-based selection at given variance threshold
    disp(['Computing ICA with the number of components which explained ' num2str(100*total_varexp) '% of the variance in PCA'])
    [~, ~, latent] = pca(data_demeaned'); % obs=time, vars=voxels
    explained = cumsum(latent) / sum(latent);
    m = find(explained >= total_varexp, 1, 'first');
    if isempty(m), m = numel(latent); end
end

%% --------------------------- Run ICA (JADE) ------------------------------

iVecs = jadeR(data_demeaned, m); % IC × voxel unmixing matrix

%% -------- Compute IC time series for each condition ----------------------

if ndims(data_temp) == 2 % no single participants nor conditions
    conds = 1;
else % several conditions
    conds = size(data_temp,3);
end
if ndims(data_temp) == 4 %single participants
    icScores = zeros(m, size(data_temp,2), conds, size(data_temp,4));
    for parti = 1:size(data_temp,4) %over participants
        for ii = 1:conds %over conditions
            icScores(:,:,ii,parti) = iVecs * data_temp(:,:,ii,parti);
        end
    end
else % provided several conditions but average across participants
    icScores = zeros(m, size(data_temp,2), conds);
    for ii = 1:conds
        icScores(:,:,ii) = iVecs * data_temp(:,:,ii);
    end
end

%% -------- Compute spatial activation patterns ----------------------------

C = cov(data'); % voxel × voxel covariance (auto-centers)
ActivationPatterns = zeros(size(data,1), m);
for c = 1:m
    ActivationPatterns(:,c) = (iVecs(c,:) * C)';
end

%% --------------------------- Pack output ---------------------------------

BROADNESS_ICA.TimeSeries_BrainNetworks         = icScores;
BROADNESS_ICA.ActivationPatterns_BrainNetworks = ActivationPatterns;
BROADNESS_ICA.numICAcomps                      = m;
BROADNESS_ICA.Time                             = time;
BROADNESS_ICA.OriginalData                     = data_temp;

%% ---------------- Helper: parse name-value pairs -------------------------

function opts = parse_name_value_pairs(opts, varargin)
    if mod(length(varargin), 2) ~= 0
        error('Arguments must be given as name-value pairs.');
    end
    for i = 1:2:length(varargin)
        name = lower(varargin{i});
        if isfield(opts, name)
            opts.(name) = varargin{i+1};
        else
            error(['Unrecognized argument: ', name]);
        end
    end
end

end

