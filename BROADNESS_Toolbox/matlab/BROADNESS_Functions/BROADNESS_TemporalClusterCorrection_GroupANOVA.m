function [significant, corrected_p, critical_p, CLUSTER_STATS] = ...
    BROADNESS_TemporalClusterCorrection_GroupANOVA(network_time_series, group_labels, alpha, varargin)
% BROADNESS_TemporalClusterCorrection_GroupANOVA
% Overall independent-group Welch ANOVA with temporal cluster-mass WILD BOOTSTRAP.
%
% INPUT
%   network_time_series : time x selected networks x conditions x participants
%   group_labels        : one integer label (1,2,...) per participant
%   alpha               : corrected cluster significance level (e.g. 0.05)
% OPTIONS
%   'permutations'      : bootstrap replicates (default 5000)
%   'clusteralpha'      : one-sided F cluster-forming level (default 0.05)
%   'randomseed'        : reproducible random seed (default 1)
%   'time'              : optional increasing time vector (seconds)
%
% OUTPUT DIMENSIONS: time x selected networks x conditions
%   significant         : samples belonging to significant clusters
%   corrected_p         : corrected cluster p repeated on observed clusters;
%                         NaN outside observed clusters (NOT pointwise p)
%   critical_p          : alpha (NOT an FDR critical p)
%   CLUSTER_STATS       : Welch F maps, chronological significant-cluster
%                         tables, and bootstrap maximum cluster masses.
%
% CORRECTION FAMILY: all time samples, SEPARATELY for each network and condition.
% NO correction across networks or conditions.
%
% Each group is centered at its observed mean. The entire residual time
% series of each participant is multiplied by one random +/-1 sign for
% each bootstrap replicate. The Welch F map is recomputed. This is an
% APPROXIMATE heteroscedasticity-aware wild bootstrap, NOT an exact
% permutation test. It requires independent participants, and its
% false-positive calibration must be validated for small/unbalanced groups.
% A significant cluster supports cluster-level inference, not individual
% sample significance or precise onset/offset inference.
% Requires Statistics and Machine Learning Toolbox (finv).

opts = struct('permutations',5000,'clusteralpha',0.05, ...
    'randomseed',1,'time',[]);
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
    error('network_time_series must be real, numeric, and finite.');
end
[nTime,nNetworks,nConditions,nParticipants] = size(network_time_series);
if nTime < 2 || nNetworks < 1 || nConditions < 1
    error('Need at least 2 time points, 1 network and 1 condition.');
end
if ~isnumeric(group_labels) || numel(group_labels) ~= nParticipants || ...
        any(~isfinite(group_labels(:))) || any(group_labels(:) < 1) || ...
        any(group_labels(:) ~= floor(group_labels(:)))
    error('group_labels must contain one positive integer per participant.');
end
group_labels = group_labels(:);
group_ids = unique(group_labels);
nGroups = numel(group_ids);
if nGroups < 3
    error('Overall group Welch ANOVA requires at least 3 groups.');
end
group_indices = cell(nGroups,1);
for g = 1:nGroups
    group_indices{g} = find(group_labels == group_ids(g));
    if numel(group_indices{g}) < 2
        error('Every group must contain at least 2 participants.');
    end
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
    error('randomseed must be an integer.');
end
if ~isempty(opts.time)
    if ~isnumeric(opts.time) || ~isreal(opts.time) || ...
            ~isvector(opts.time) || numel(opts.time) ~= nTime || ...
            any(~isfinite(opts.time(:))) || any(diff(opts.time(:)) <= 0)
        error('time must be a finite, strictly increasing vector of nTime values.');
    end
    opts.time = opts.time(:);
end

significant = false(nTime,nNetworks,nConditions);
corrected_p = nan(nTime,nNetworks,nConditions);
critical_p = alpha;
f_statistics = zeros(nTime,nNetworks,nConditions);
cluster_labels = zeros(nTime,nNetworks,nConditions);
cluster_tables = cell(nNetworks,nConditions);
null_maximum_mass = zeros(opts.permutations,nNetworks,nConditions);
cluster_threshold = nan(nTime,nNetworks,nConditions);

previous_rng = rng;
cleanup_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(opts.randomseed,'twister');

for net = 1:nNetworks
    for condition = 1:nConditions
        group_data = cell(nGroups,1);
        group_residuals = cell(nGroups,1);
        for g = 1:nGroups
            group_data{g} = reshape(network_time_series(:,net,condition,group_indices{g}),nTime,[]);
            group_residuals{g} = bsxfun(@minus,group_data{g},mean(group_data{g},2));
        end
        [observed_f,observed_df2] = welch_f_map(group_data);
        observed_threshold = finv(1-opts.clusteralpha,nGroups-1,observed_df2);
        f_statistics(:,net,condition) = observed_f;
        cluster_threshold(:,net,condition) = observed_threshold;
        [labels,masses] = temporal_f_clusters(observed_f,observed_threshold);
        cluster_labels(:,net,condition) = labels;

        for perm = 1:opts.permutations
            boot_data = cell(nGroups,1);
            for g = 1:nGroups
                signs = 2*(rand(1,size(group_residuals{g},2)) > 0.5)-1;
                boot_data{g} = bsxfun(@times,group_residuals{g},signs);
            end
            [boot_f,boot_df2] = welch_f_map(boot_data);
            boot_threshold = finv(1-opts.clusteralpha,nGroups-1,boot_df2);
            [~,boot_masses] = temporal_f_clusters(boot_f,boot_threshold);
            if ~isempty(boot_masses)
                null_maximum_mass(perm,net,condition) = max(boot_masses);
            end
        end

        cluster_p = nan(numel(masses),1);
        starts = zeros(numel(masses),1);
        ends = zeros(numel(masses),1);
        mean_stat = nan(numel(masses),1);
        peak_stat = nan(numel(masses),1);
        for c = 1:numel(masses)
            members = find(labels == c);
            starts(c) = members(1);
            ends(c) = members(end);
            mean_stat(c) = mean(observed_f(members));
            peak_stat(c) = max(observed_f(members));
            cluster_p(c) = (1+sum(null_maximum_mass(:,net,condition) >= masses(c))) ...
                / (opts.permutations+1);
            corrected_p(members,net,condition) = cluster_p(c);
            significant(members,net,condition) = cluster_p(c) <= alpha;
        end
        keep = find(cluster_p <= alpha);
        [~,order] = sort(starts(keep));
        keep = keep(order);
        cluster_number = (1:numel(keep))';
        if isempty(opts.time)
            cluster_tables{net,condition} = table(cluster_number, ...
                starts(keep),ends(keep),mean_stat(keep),peak_stat(keep), ...
                masses(keep),cluster_p(keep), ...
                'VariableNames',{'Cluster','StartSample','EndSample', ...
                'MeanStatistic','PeakStatistic','Mass','CorrectedP'});
        else
            cluster_tables{net,condition} = table(cluster_number, ...
                opts.time(starts(keep)),opts.time(ends(keep)), ...
                mean_stat(keep),peak_stat(keep),masses(keep),cluster_p(keep), ...
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
CLUSTER_STATS.Settings.GroupLabels = group_labels;
CLUSTER_STATS.Settings.GroupIDs = group_ids;
CLUSTER_STATS.Settings.CorrectionFamily = ...
    'All time points, separately per network and condition';
CLUSTER_STATS.Settings.Resampling = ...
    'Within-group centered-residual Rademacher wild bootstrap';
CLUSTER_STATS.Settings.InferenceCaveat = ...
    'Approximate cluster-level inference; validate with null simulations, especially for small or unequal groups';
CLUSTER_STATS.FStatistics = f_statistics;
CLUSTER_STATS.ClusterFormingThreshold = cluster_threshold;
CLUSTER_STATS.ClusterLabels = cluster_labels;
CLUSTER_STATS.Clusters = cluster_tables;
CLUSTER_STATS.NullMaximumClusterMass = null_maximum_mass;
CLUSTER_STATS.ClusterCorrectedP = corrected_p;
CLUSTER_STATS.Significant = significant;
end

function [f,df2] = welch_f_map(group_data)
    k = numel(group_data);
    nTime = size(group_data{1},1);
    n = zeros(1,k);
    m = zeros(nTime,k);
    v = zeros(nTime,k);
    for g = 1:k
        x = group_data{g};
        n(g) = size(x,2);
        m(:,g) = mean(x,2);
        v(:,g) = var(x,0,2);
    end
    % Welch's ANOVA is undefined for zero within-group variance.
    % Such samples are excluded from cluster formation.
    valid = all(isfinite(v) & v > 0,2);
    f = zeros(nTime,1);
    df2 = nan(nTime,1);
    if ~any(valid)
        return
    end
    w = bsxfun(@rdivide,n,v(valid,:));
    wsum = sum(w,2);
    weighted_mean = sum(w.*m(valid,:),2)./wsum;
    correction_sum = sum(bsxfun(@rdivide,(1-bsxfun(@rdivide,w,wsum)).^2,n-1),2);
    df1 = k-1;
    df2(valid) = (k^2-1)./(3*correction_sum);
    f(valid) = (sum(w.*bsxfun(@minus,m(valid,:),weighted_mean).^2,2)/df1) ...
        ./ (1+2*(k-2)*correction_sum/(k^2-1));
    f(~isfinite(f) | f < 0) = 0;
end

function [labels,masses] = temporal_f_clusters(f,threshold)
    active = isfinite(threshold(:)) & f(:) >= threshold(:);
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
