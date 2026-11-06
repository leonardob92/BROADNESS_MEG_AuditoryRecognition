function [significant, corrected_p, critical_p, CLUSTER_STATS] = BROADNESS_TemporalClusterCorrection_ANOVA(network_time_series, alpha, varargin)
% BROADNESS_TemporalClusterCorrection_ANOVA
% Temporal cluster-mass permutation test for the omnibus within-participant
% condition effect, separately for each selected brain network.
%
% INPUT
%   network_time_series : time x selected networks x conditions x participants
%   alpha               : cluster-level significance threshold (e.g. 0.05)
%
% OPTIONAL NAME-VALUE INPUTS
%   'permutations'       : random permutations (default 5000)
%   'clusteralpha'       : one-sided F cluster-forming p (default 0.05)
%   'randomseed'         : reproducible random seed (default 1)
%   'time'              : optional time vector in seconds (default [])
%
% OUTPUT (time x selected networks)
%   significant         : true at samples in significant temporal clusters
%   corrected_p         : CLUSTER-LEVEL corrected p repeated across samples
%                         of each observed suprathreshold cluster; NaN elsewhere
%   critical_p          : alpha (not an FDR critical p-value)
%   CLUSTER_STATS       : observed F maps, cluster tables, null maxima, settings
%                         Tables contain only significant clusters, sorted
%                         chronologically and numbered consecutively.
%                         Time boundaries are used if 'time' is supplied;
%                         otherwise sample indices are reported.
%                         MeanStatistic and PeakStatistic are F-statistics.
%
% PERMUTATION PROCEDURE
%   Within each participant, permute complete condition time series (one
%   permutation of condition labels per participant, shared across time).
%   Recompute the repeated-measures F at each time point. Record the maximum
%   positive-F cluster mass across ALL time points in each permutation.
%   Each network has its OWN permutation null distribution. No multiplicity
%   correction is applied across networks.
%
% ASSUMPTIONS AND INTERPRETATION
%   The within-participant condition labels must be exchangeable under the
%   omnibus null. This may not hold for arbitrary condition covariance
%   structures or designs with ordering/carryover effects. The statistic
%   is the classical repeated-measures F, calibrated by permutation (not a
%   Greenhouse-Geisser-adjusted parametric ANOVA). Significant clusters
%   support CLUSTER-level inference, not significance at each sample or
%   precise onset/offset times. This is not a correction for pairwise tests.
%
% Requires Statistics and Machine Learning Toolbox (finv).

opts.permutations = 5000;
opts.clusteralpha = 0.05;
opts.randomseed = 1;
opts.time = [];
if mod(numel(varargin),2) ~= 0
    error('Optional arguments must be name-value pairs.');
end
for k = 1:2:numel(varargin)
    key = lower(char(varargin{k}));
    if ~isfield(opts,key)
        error('Unknown option: %s',key);
    end
    opts.(key) = varargin{k+1};
end

if ~isnumeric(network_time_series) || ~isreal(network_time_series) || ...
        any(~isfinite(network_time_series(:)))
    error('network_time_series must be real, numeric and finite.');
end
[nTime,nNetworks,nConditions,nParticipants] = size(network_time_series);
if ~isempty(opts.time)
    if ~isnumeric(opts.time) || ~isreal(opts.time) || ...
            ~isvector(opts.time) || numel(opts.time) ~= nTime || ...
            any(~isfinite(opts.time(:)))
        error('time must be a finite numeric vector with one value per time point.');
    end
    opts.time = opts.time(:);
end
if nTime < 2 || nNetworks < 1 || nConditions < 2 || nParticipants < 2
    error('Require >=2 time points, >=1 network, >=2 conditions and >=2 participants.');
end
if ~isnumeric(alpha) || ~isscalar(alpha) || ~isfinite(alpha) || ...
        alpha <= 0 || alpha >= 1 || ...
        ~isnumeric(opts.clusteralpha) || ~isscalar(opts.clusteralpha) || ...
        ~isfinite(opts.clusteralpha) || opts.clusteralpha <= 0 || opts.clusteralpha >= 1
    error('alpha and clusteralpha must be numbers strictly between 0 and 1.');
end
if ~isnumeric(opts.permutations) || ~isscalar(opts.permutations) || ...
        ~isfinite(opts.permutations) || opts.permutations < 1 || ...
        opts.permutations ~= floor(opts.permutations)
    error('permutations must be a positive integer.');
end
if ~isnumeric(opts.randomseed) || ~isscalar(opts.randomseed) || ...
        ~isfinite(opts.randomseed) || opts.randomseed ~= floor(opts.randomseed)
    error('randomseed must be a finite integer.');
end

% Standard repeated-measures ANOVA degrees of freedom.
dfCondition = nConditions-1;
dfError = (nParticipants-1)*dfCondition;
fThreshold = finv(1-opts.clusteralpha,dfCondition,dfError);

significant = false(nTime,nNetworks);
corrected_p = nan(nTime,nNetworks);
critical_p = alpha;
f_statistics = zeros(nTime,nNetworks);
cluster_labels = zeros(nTime,nNetworks);
cluster_tables = cell(nNetworks,1);
null_maximum_mass = zeros(opts.permutations,nNetworks);

previous_rng = rng;
cleanup_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(opts.randomseed,'twister');

for net = 1:nNetworks
    % Explicit reshape preserves the condition and participant dimensions.
    x = reshape(network_time_series(:,net,:,:),nTime,nConditions,nParticipants);
    observedF = repeated_f(x);
    f_statistics(:,net) = observedF;
    [labels,masses] = positive_clusters(observedF,fThreshold);
    cluster_labels(:,net) = labels;

    for perm = 1:opts.permutations
        permuted = zeros(size(x));
        for participant = 1:nParticipants
            order = randperm(nConditions);
            permuted(:,:,participant) = x(:,order,participant);
        end
        permutedF = repeated_f(permuted);
        [~,permutedMasses] = positive_clusters(permutedF,fThreshold);
        if ~isempty(permutedMasses)
            null_maximum_mass(perm,net) = max(permutedMasses);
        end
    end

    start_sample = zeros(numel(masses),1);
    end_sample = zeros(numel(masses),1);
    mean_statistic = nan(numel(masses),1);
    peak_statistic = nan(numel(masses),1);
    cluster_p = nan(numel(masses),1);
    for c = 1:numel(masses)
        members = find(labels == c);
        start_sample(c) = members(1);
        end_sample(c) = members(end);
        cluster_f = observedF(members);
        mean_statistic(c) = mean(cluster_f);
        peak_statistic(c) = max(cluster_f);
        cluster_p(c) = (1+sum(null_maximum_mass(:,net) >= masses(c))) / ...
            (opts.permutations+1);
        corrected_p(members,net) = cluster_p(c);
        significant(members,net) = cluster_p(c) <= alpha;
    end
    % Report only significant clusters, ordered by their onset.
    % Original cluster IDs remain unchanged in ClusterLabels.
    keep_indices = find(cluster_p <= alpha);
    if isempty(opts.time)
        [~,order] = sort(start_sample(keep_indices));
    else
        [~,order] = sort(opts.time(start_sample(keep_indices)));
    end
    keep_indices = keep_indices(order);
    display_id = (1:numel(keep_indices))';
    if isempty(opts.time)
        cluster_tables{net} = table(display_id, ...
            start_sample(keep_indices),end_sample(keep_indices), ...
            mean_statistic(keep_indices),peak_statistic(keep_indices), ...
            masses(keep_indices),cluster_p(keep_indices), ...
            'VariableNames',{'Cluster','StartSample','EndSample', ...
            'MeanStatistic','PeakStatistic','Mass','CorrectedP'});
    else
        start_time = opts.time(start_sample);
        end_time = opts.time(end_sample);
        cluster_tables{net} = table(display_id, ...
            start_time(keep_indices),end_time(keep_indices), ...
            mean_statistic(keep_indices),peak_statistic(keep_indices), ...
            masses(keep_indices),cluster_p(keep_indices), ...
            'VariableNames',{'Cluster','StartTime','EndTime', ...
            'MeanStatistic','PeakStatistic','Mass','CorrectedP'});
    end
end

CLUSTER_STATS.Settings.Alpha = alpha;
CLUSTER_STATS.Settings.ClusterFormingAlpha = opts.clusteralpha;
CLUSTER_STATS.Settings.Permutations = opts.permutations;
CLUSTER_STATS.Settings.RandomSeed = opts.randomseed;
CLUSTER_STATS.Settings.Time = opts.time;
CLUSTER_STATS.Settings.CorrectionFamily = ...
    'All time points, separately per network (omnibus condition effect)';
CLUSTER_STATS.Settings.Resampling = ...
    'Complete condition-label permutations within participants';
CLUSTER_STATS.Settings.InferenceCaveat = ...
    'Condition exchangeability required; cluster-level inference only; no across-network correction';
CLUSTER_STATS.FStatistics = f_statistics;
CLUSTER_STATS.ClusterFormingThreshold = fThreshold;
CLUSTER_STATS.DegreesFreedom = [dfCondition dfError];
CLUSTER_STATS.ClusterLabels = cluster_labels;
CLUSTER_STATS.Clusters = cluster_tables;
CLUSTER_STATS.NullMaximumClusterMass = null_maximum_mass;
CLUSTER_STATS.ClusterCorrectedP = corrected_p;
CLUSTER_STATS.Significant = significant;
end

function f = repeated_f(x)
% x: time x conditions x participants. Classical one-factor RM F.
    [~,nConditions,nParticipants] = size(x);
    grand = mean(mean(x,3),2);
    conditionMeans = mean(x,3);
    participantMeans = mean(x,2);
    ssCondition = nParticipants * sum((conditionMeans-grand).^2,2);
    ssParticipant = nConditions * sum((participantMeans-grand).^2,3);
    ssTotal = sum(sum((x-grand).^2,3),2);
    ssError = max(0,ssTotal-ssCondition-ssParticipant);
    df1 = nConditions-1;
    df2 = (nParticipants-1)*df1;
    f = (ssCondition/df1) ./ (ssError/df2);
    % If both effect and error are zero, there is no observed effect.
    f(isnan(f)) = 0;
    f = f(:);
end

function [labels,masses] = positive_clusters(f,threshold)
    active = f(:) >= threshold;
    changes = diff([false;active;false]);
    starts = find(changes == 1);
    ends = find(changes == -1)-1;
    labels = zeros(numel(f),1);
    masses = zeros(numel(starts),1);
    for c = 1:numel(starts)
        members = starts(c):ends(c);
        labels(members) = c;
        masses(c) = sum(f(members));
    end
end
