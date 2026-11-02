function [significant, corrected_p, critical_p, CLUSTER_STATS] = BROADNESS_TemporalClusterCorrection(network_time_series, condition_pairs, alpha, varargin)
% BROADNESS_TemporalClusterCorrection
% Paired temporal cluster-mass permutation test for within-participant conditions.
%
% INPUT
%   network_time_series : time x selected networks x conditions x participants
%   condition_pairs     : Npairs x 2 matrix of condition indices
%   alpha               : corrected cluster significance level (e.g. 0.05)
%
% OPTIONAL NAME-VALUE INPUTS
%   'permutations'      : number of random permutations (default 5000)
%   'clusteralpha'      : two-sided cluster-forming p threshold (default 0.05)
%   'randomseed'        : seed for reproducibility (default 1)
%   'time'              : optional time vector in seconds (default [])
%
% OUTPUT (time x selected networks x condition pairs)
%   significant         : true at samples belonging to significant clusters
%   corrected_p         : CLUSTER-LEVEL corrected p, repeated across each
%                         observed supra-threshold cluster; NaN elsewhere.
%                         Do not interpret as pointwise adjusted p-values.
%   critical_p          : alpha, the corrected cluster significance threshold
%   CLUSTER_STATS       : details, including t maps, cluster tables, and null maxima
%                         NullMaximumClusterMass: permutations x networks x pairs
%                         Clusters contains only significant clusters.
%                         Time boundaries are used if 'time' is supplied;
%                         otherwise sample indices are reported.
%                         Reported clusters are chronological and renumbered
%                         1,2,...; ClusterLabels retains original cluster IDs.
%
% Each network AND each selected condition pair has its own null maximum
% cluster-mass distribution across ALL time samples. Positive and negative
% clusters are kept separate, but share a two-sided maximum-mass null.
%
% Null resampling: randomly flip the sign of each participant's COMPLETE
% paired-difference time series (+1 or -1). This preserves the temporal
% dependence within participants and assumes sign-symmetric differences
% under the null. Grouped correction across condition pairs or networks
% is NOT performed; plan contrasts before examining the results.
%
% Cluster inference concerns whole clusters, not individual sample onsets.
% Requires Statistics and Machine Learning Toolbox (tinv).

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
if nTime < 2 || nNetworks < 1 || nConditions < 2 || nParticipants < 2
    error('Require >=2 time points, >=1 network, >=2 conditions and >=2 participants.');
end
if ~isempty(opts.time)
    if ~isnumeric(opts.time) || ~isreal(opts.time) || ...
            ~isvector(opts.time) || numel(opts.time) ~= nTime || ...
            any(~isfinite(opts.time(:)))
        error('time must be a finite numeric vector with one value per time point.');
    end
    opts.time = opts.time(:);
end
if ~isnumeric(condition_pairs) || size(condition_pairs,2) ~= 2 || ...
        isempty(condition_pairs) || any(~isfinite(condition_pairs(:))) || ...
        any(condition_pairs(:) ~= floor(condition_pairs(:))) || ...
        any(condition_pairs(:) < 1) || any(condition_pairs(:) > nConditions) || ...
        any(condition_pairs(:,1) == condition_pairs(:,2)) || ...
        size(unique(sort(condition_pairs,2),'rows'),1) ~= size(condition_pairs,1)
    error('condition_pairs must contain distinct, valid, nonduplicated pairs.');
end
if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 0 || alpha >= 1 || ...
        ~isscalar(opts.clusteralpha) || ~isfinite(opts.clusteralpha) || ...
        opts.clusteralpha <= 0 || opts.clusteralpha >= 1
    error('alpha and clusteralpha must be between 0 and 1.');
end
if ~isscalar(opts.permutations) || ~isfinite(opts.permutations) || ...
        opts.permutations < 1 || opts.permutations ~= floor(opts.permutations)
    error('permutations must be a positive integer.');
end
if ~isscalar(opts.randomseed) || ~isfinite(opts.randomseed) || ...
        opts.randomseed ~= floor(opts.randomseed)
    error('randomseed must be a finite integer.');
end

nPairs = size(condition_pairs,1);
significant = false(nTime,nNetworks,nPairs);
corrected_p = nan(nTime,nNetworks,nPairs);
critical_p = alpha;
t_statistics = zeros(nTime,nNetworks,nPairs);
cluster_labels = zeros(nTime,nNetworks,nPairs);
cluster_tables = cell(nNetworks,nPairs);
null_maximum_mass = zeros(opts.permutations,nNetworks,nPairs);
threshold = tinv(1-opts.clusteralpha/2,nParticipants-1);

previous_rng = rng;
cleanup_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(opts.randomseed,'twister');

for net = 1:nNetworks
    % Explicit reshape keeps the condition dimension when nConditions = 2.
    x = reshape(network_time_series(:,net,:,:), ...
        nTime,nConditions,nParticipants);
    observed_masses = cell(nPairs,1);
    observed_labels = cell(nPairs,1);
    for pair = 1:nPairs
        d = reshape(x(:,condition_pairs(pair,1),:) - ...
            x(:,condition_pairs(pair,2),:),nTime,nParticipants);
        t = paired_t(d);
        t_statistics(:,net,pair) = t;
        [labels,masses] = temporal_clusters(t,threshold);
        observed_labels{pair} = labels;
        observed_masses{pair} = masses;
        cluster_labels(:,net,pair) = labels;
    end

    % Each contrast has an independent sign-flip null distribution.
    for pair = 1:nPairs
        d = reshape(x(:,condition_pairs(pair,1),:) - ...
            x(:,condition_pairs(pair,2),:),nTime,nParticipants);
        for perm = 1:opts.permutations
            signs = 2*(rand(1,nParticipants) > 0.5)-1;
            permuted_d = bsxfun(@times,d,signs);
            permuted_t = paired_t(permuted_d);
            [~,masses] = temporal_clusters(permuted_t,threshold);
            if ~isempty(masses)
                null_maximum_mass(perm,net,pair) = max(masses);
            end
        end
    end

    for pair = 1:nPairs
        labels = observed_labels{pair};
        masses = observed_masses{pair};
        cluster_id = (1:numel(masses))';
        cluster_p = nan(numel(masses),1);
        start_sample = zeros(numel(masses),1);
        end_sample = zeros(numel(masses),1);
        mean_statistic = nan(numel(masses),1);
        peak_statistic = nan(numel(masses),1);
        for c = 1:numel(masses)
            members = find(labels == c);
            start_sample(c) = members(1);
            end_sample(c) = members(end);
            cluster_t = t_statistics(members,net,pair);
            mean_statistic(c) = mean(cluster_t);
            [~,peak_index] = max(abs(cluster_t));
            peak_statistic(c) = cluster_t(peak_index);
            cluster_p(c) = (1+sum(null_maximum_mass(:,net,pair) >= masses(c))) / ...
                (opts.permutations+1);
            corrected_p(members,net,pair) = cluster_p(c);
            significant(members,net,pair) = cluster_p(c) <= alpha;
        end
        % Report only significant clusters, ordered by their onset.
        % The original cluster IDs remain unchanged in ClusterLabels.
        keep_indices = find(cluster_p <= alpha);
        if isempty(opts.time)
            [~,order] = sort(start_sample(keep_indices));
        else
            [~,order] = sort(opts.time(start_sample(keep_indices)));
        end
        keep_indices = keep_indices(order);
        display_id = (1:numel(keep_indices))';
        if isempty(opts.time)
            cluster_tables{net,pair} = table(display_id, ...
                start_sample(keep_indices),end_sample(keep_indices), ...
                mean_statistic(keep_indices),peak_statistic(keep_indices), ...
                masses(keep_indices),cluster_p(keep_indices), ...
                'VariableNames',{'Cluster','StartSample','EndSample', ...
                'MeanStatistic','PeakStatistic','Mass','CorrectedP'});
        else
            start_time = opts.time(start_sample);
            end_time = opts.time(end_sample);
            cluster_tables{net,pair} = table(display_id, ...
                start_time(keep_indices),end_time(keep_indices), ...
                mean_statistic(keep_indices),peak_statistic(keep_indices), ...
                masses(keep_indices),cluster_p(keep_indices), ...
                'VariableNames',{'Cluster','StartTime','EndTime', ...
                'MeanStatistic','PeakStatistic','Mass','CorrectedP'});
        end
    end
end

CLUSTER_STATS.Settings.Alpha = alpha;
CLUSTER_STATS.Settings.ClusterFormingAlpha = opts.clusteralpha;
CLUSTER_STATS.Settings.Permutations = opts.permutations;
CLUSTER_STATS.Settings.RandomSeed = opts.randomseed;
CLUSTER_STATS.Settings.Time = opts.time;
CLUSTER_STATS.Settings.ConditionPairs = condition_pairs;
CLUSTER_STATS.Settings.CorrectionFamily = ...
    'All time points, separately per network and condition pair';
CLUSTER_STATS.Settings.Resampling = 'Paired-difference participant sign flips';
CLUSTER_STATS.Settings.InferenceCaveat = ...
    'Sign-symmetric paired differences required; no correction across pairs or networks';
CLUSTER_STATS.TStatistics = t_statistics;
CLUSTER_STATS.ClusterFormingThreshold = threshold;
CLUSTER_STATS.ClusterLabels = cluster_labels;
CLUSTER_STATS.Clusters = cluster_tables;
CLUSTER_STATS.NullMaximumClusterMass = null_maximum_mass;
CLUSTER_STATS.ClusterCorrectedP = corrected_p;
CLUSTER_STATS.Significant = significant;
end

function t = paired_t(d)
    n = size(d,2);
    avg = mean(d,2);
    sd = std(d,0,2);
    t = avg ./ (sd ./ sqrt(n));
    % Constant zero differences have no evidence of an effect. Nonzero
    % constant differences produce +/-Inf, which form suprathreshold clusters.
    t(isnan(t)) = 0;
end

function [labels,masses] = temporal_clusters(t,threshold)
    labels = zeros(numel(t),1);
    masses = zeros(0,1);
    for direction = [1 -1]
        active = (direction*t(:)) >= threshold;
        changes = diff([false;active;false]);
        starts = find(changes == 1);
        ends = find(changes == -1)-1;
        for c = 1:numel(starts)
            members = starts(c):ends(c);
            cluster = numel(masses)+1;
            labels(members) = cluster;
            masses(cluster,1) = sum(abs(t(members))); %#ok<AGROW>
        end
    end
end
