function [BROADNESS] = BROADNESS_NetworkEstimation(data, time, varargin)

% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  NETWORK ESTIMATION
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
%  This function takes a multivariate dataset (channels or brain voxels × time matrix [x experimental conditions x participants])
%  and returns the estimated broadband brain networks.
%  BROADNESS is particularly suitable for event-related design and for
%  estimating brain networks in the context of event-related fiels (ERF) or
%  potentials (ERP).
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - data   :  DATA MATRIX, either:
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
%  - time   : vector with time in seconds
%
%
%  The following parameters can be specified using name-value pairs:
%
%  - 'time_window'      : Segment of the data to be used in the analysis (vector in seconds, e.g. [0.350 1.750]).
%                         Default: full data duration.
%
%  - 'permutations_num' : Number of permutations for Monte-Carlo simulations (MCS)
%                         (e.g. 100) to establish how many brain networks (PCs) are significant.
%                         Default: 0 (no MCS computed)
%
%  - 'randomization'    : Strategy for randomisation (relevant only if permutations_num is provided as a positive integer) 
%                            1 = randomising only time-points (independently for each time series)
%                            2 = randomising only space (independently for each time-point) 
%                            3 = randomising both time and space
%                         Default: 1
%
%  - 'sign_eigenvect'   : Method for normalizing eigenvectors sign (SUGGESTED EITHER 'max_abs' or 'average'):
%                            -'occurrences' = using mean of the negative/positive values occurrences
%                            -'max_abs' = on the basis of the sign of the maximum value in absolute terms
%                            -'average' = on the basis of the sign of the average of weights
%                         Default: max_abs
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - BROADNESS   : Structure containing the following information for each brain network (PC):
%                   1. spatial activation patterns (eigenvectors)
%                   2. variance explained (eigenvalues) - original data
%                   3. brain network time series
%                   4. number of significant brain networks after MCS
%                   5. variance explained (eigenvalues) - permuted data (only if MCS was computed) 
%                   6. time (carried out for future plotting purposes) 
%                   7. data (carried out for future plotting purposes)
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






%% handling optional arguments

disp('Checking inputs')

% Defining default values
opts = struct('time_window', [], 'permutations_num', 0, 'randomization', 1, 'sign_eigenvect', 'max_abs');

% Parsing name-value pair arguments
opts = parse_name_value_pairs(opts, varargin{:});

% Assigning values
permutations_num = opts.permutations_num;
randomization = opts.randomization;
sign_eigenvect = opts.sign_eigenvect;

% Checking inputs
if randomization < 1 || randomization > 3
    error('Randomization must be between 1 and 3')
end
if permutations_num < 0
    warning('Your permutations_num is negative.. and thus it does not compute anything; are you sure this is what you want?')
end
sz = size(data);
non_singleton_dims = sum(sz > 1); %trick to get if the matrix is a vector
if non_singleton_dims < 2 || non_singleton_dims > 4
    error('The data matrix must be either 2D (e.g. brain voxels x time), 3D (brain voxels x time x experimental conditions), or 4D (brain voxels x time x experimental conditions x participants)')
end
if size(data,2) ~= length(time)
    error('The number of time-points in the data matrix (2nd dimension) must be equal to the time-points in the time vector')
end
    
% If duration is empty, use full data duration
if isempty(opts.time_window)
    disp('Duration not provided.. using full data duration');
    idx_start = 1;
    idx_end = size(data,2);
else
    % Find indices of closest time points
    [~, idx_start] = min(abs(time - opts.time_window(1)));
    [~, idx_end]   = min(abs(time - opts.time_window(2)));
    if idx_end <= idx_start
        error('The end of your time-window cannot be smaller than the beginning')
    end
end


%% PCA computation

disp('Computing PCA')

if non_singleton_dims == 4 %if data is provided for several independent participants
    data_temp = data; %store the data
    data = mean(data,4); %average it across participants
end

data = data(:,idx_start:idx_end,:); %extracting time-window of interest
time = time(idx_start:idx_end); %same for time vector in seconds
averaged_data = mean(data(:,:,:),3); %average across 3rd dimension (e.g. experimental conditions)

if exist('pca', 'file') == 2
    % The PCA function exists (from Statistics Toolbox)
    disp('Using built-in PCA function...');
    [activation_patterns,~,~,~,variance] = pca(averaged_data'); %actual computation of PCA
else
    % Fallback if PCA function is not available
    disp('pca() not found — using custom PCA implementation...');
    data_demeaned = bsxfun(@minus,averaged_data,mean(averaged_data)); %demeaning of the data
    data_covariance = cov(data_demeaned'); %covariance matrix (note that data_demeaned is transposed, so time-points x brain sources)
    [activation_patterns,eigenvalues] = eig(data_covariance); %eigenvector solution
    [eigenvalues,sidx]  = sort( diag(eigenvalues),'descend' ); % the first output returns sorted evals extracted from diagonal
    activation_patterns = activation_patterns(:,sidx);          % sort eigenvectors
    variance = eigenvalues.*100./sum(eigenvalues); % normalize eigenvalues to percent variance explained
end


%% PCA on randomized data

if permutations_num > 0 %if MCS was requested
    
    disp('Computing PCA on randomized data')
    
    max_variance_permutations = zeros(permutations_num,1); %preallocating vector for variance explained by 1st PC for each permutation
    variance_randomized_perms = []; %variance of the permuted data to be stored for future plotting purposes
    for permi = 1:permutations_num %over permutations
        if randomization == 1     %randomizing only time
            data_reshaped = zeros(size(averaged_data,1),size(averaged_data,2)); %preallocating matrix for randomized data
            for sourci = 1:size(averaged_data,1) %over brain sources
                idx_dummy = randperm(size(averaged_data,2)); %creating a permuted array of the indices of the temporal dimension
                data_reshaped(sourci,:) = averaged_data(sourci,idx_dummy); %actual reshaped data
            end
        elseif randomization == 2 %randomizing only space
            data_reshaped = zeros(size(averaged_data,1),size(averaged_data,2)); %preallocating matrix for randomized data
            for sourci = 1:size(averaged_data,2) %over time-points
                idx_dummy = randperm(size(averaged_data,1)); %creating a permuted array of the indices of the spatial dimension
                data_reshaped(:,sourci) = averaged_data(idx_dummy,sourci); %actual reshaped data
            end
        elseif randomization == 3 %randomizing both space and time
            idx_dummy = randperm(size(averaged_data,1)*size(averaged_data,2)); %creating a permuted array from 1 to size of original data 
            data_reshaped_dummy = zeros(size(averaged_data,1)*size(averaged_data,2),1); %initialising new temporary vector
            data_reshaped_dummy(1:size(averaged_data,1)*size(averaged_data,2)) = averaged_data(idx_dummy); %taking elements in matrix data with index idx_dummy
            data_reshaped = reshape(data_reshaped_dummy,[size(averaged_data,1),size(averaged_data,2)]); %reshaping the vector into a matrix shaped as the original data
        end
        if exist('pca', 'file') == 2
            % The PCA function exists (from Statistics Toolbox)
            [~,~,~,~,variance_randomized] = pca(data_reshaped'); %PCA on randomized data
        else
            % Fallback if PCA function is not available
            data_demeaned_r = bsxfun(@minus,data_reshaped,mean(data_reshaped)); %demeaning of the data
            data_covariance_r = cov(data_demeaned_r'); %covariance matrix (note that data_demeaned is transposed, so time-points x brain sources)
            [~,eigenvalues_r] = eig(data_covariance_r); %eigenvector solution
            [eigenvalues_r,sidx_r]  = sort( diag(eigenvalues_r),'descend' ); % the first output returns sorted evals extracted from diagonal
            variance_randomized = eigenvalues_r.*100./sum(eigenvalues_r); % normalize eigenvalues to percent variance explained
        end
        variance_randomized_perms = cat(2,variance_randomized_perms,variance_randomized); %storing variance of the permuted data to be stored for future plotting purposes
        max_variance_permutations(permi,1) = max(max(variance_randomized)); %storing maximum variance (max eigenvalue) occurring for PCs in randomized data
        disp(['Permutation number ' num2str(permi) ' / ' num2str(permutations_num)])
    end
    variance_randomized_perms = mean(variance_randomized_perms,2); %average across permutations
    max_variance = max(max_variance_permutations); %maximum variance over permutations
    PCs = find(variance>max_variance); %significant PCs obtained by getting PCs with variance (eigenvalues) higher than the maximum variance (max eigenvalues over the permutations) obtained from randomized data
    if isempty(PCs)
        disp('There are no PCs (brain networks) that survived MCS')
    else
        disp('Percentage of variance explained by significant PCs (brain networks)')
        disp(variance(PCs))
    end
else
    PCs = 1:length(variance); %keeping all PCs
end


%% normalizing eigenvectors signs

dummy_ones = ones(size(activation_patterns,1),size(activation_patterns,2)); %vector of 1s with length of significant PCs
switch sign_eigenvect
    case 'occurrences'
        % 1) normalizing eigenvectors sign by using mean of the negative/positive values occurrences
        dummy_ones(:,mean(activation_patterns> 0) < 0.5) = -1; %assigning -1 to eigenvectors that have more positive than negative weights
    case 'max_abs'
        % 2) normalizing eigenvectors sign on the basis of the sign of the maximum value in absolute terms
        [~,mi] = max(abs(activation_patterns)); %getting maximum values indices
        ab = zeros(1,length(mi));
        for ii = 1:length(mi) %over significant PCs
            ab(1,ii) = activation_patterns(mi(ii),ii); %storing original values with signs (corresponding to maximum in absolute terms)
        end
        dummy_ones(:,sign(ab)<0) = -1;
    case 'average'
        % 3) normalizing eigenvectors sign on the basis of the sign of the average of activation patterns
        mVV = mean(activation_patterns);
        dummy_ones(:,sign(mVV)<0) = -1;
end
activation_patterns = activation_patterns .* dummy_ones; %actual normalization of the eigenvectors sign


%% preparing output structure

disp('Preparing output')

BROADNESS = [];
BROADNESS.Variance_BrainNetworks = variance; %storing variance of significant brain networks (PCs)

if permutations_num > 0 %if MCS was run
    BROADNESS.Significant_BrainNetworks = PCs;
    BROADNESS.VariancePermutations = variance_randomized_perms;
else
    BROADNESS.Significant_BrainNetworks = 'MCS has not been run';
end

BROADNESS.ActivationPatterns_BrainNetworks = activation_patterns; %weights of PCA

if isempty(PCs) %very unlikely case where no PCs of the original data explain a higher variance than the randomized data
    TimeSeries = [];
else
    if non_singleton_dims == 4 %if data is provided for several independent participants
        TimeSeries = zeros(length(time),PCs(end),size(data_temp,3),size(data_temp,4));
        for parti = 1:size(data_temp,4) %over participants
            disp(['Computing time series for participant ' num2str(parti) ' / ' num2str(size(data_temp,4))])
            for condi = 1:size(data,3) %over experimental conditions (or whatever the user has in the 3rd dimension of the data matrix)
                TimeSeries(:,:,condi,parti) = data_temp(:,1:length(time),condi,parti)' * activation_patterns(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
            end
        end
    else
        TimeSeries = zeros(length(time),PCs(end),size(data,3)); %preallocating matrix for brain networks time series
        for condi = 1:size(data,3) %over experimental conditions (or whatever the user has in the 3rd dimension of the data matrix)
            TimeSeries(:,:,condi) = data(:,1:length(time),condi)' * activation_patterns(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
        end
        
    end
end

BROADNESS.TimeSeries_BrainNetworks = TimeSeries; %storing PCA time series

BROADNESS.Time = time; %storing time
if non_singleton_dims == 4 %if data is provided for several independent participants
    BROADNESS.OriginalData = data_temp; %storing original data
else
    BROADNESS.OriginalData = data;
end


%%

%% Helper Function: Parse Name-Value Pairs

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

%%

end
