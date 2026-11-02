function [significant, corrected_p, critical_p, CLUSTER_STATS] = ...
    BROADNESS_TemporalClusterCorrection_Groups(network_time_series, group_labels, group_pairs, alpha, varargin)
% BROADNESS_TemporalClusterCorrection_Groups
% Pairwise independent-group temporal cluster-mass WILD BOOTSTRAP.
%
% INPUTS
%   network_time_series : time x selected networks x conditions x participants
%   group_labels        : one integer group label (1,2,...) per participant
%   group_pairs         : group pairs (rows), e.g. nchoosek(1:3,2)
%   alpha               : cluster-level significance threshold (e.g. 0.05)
%
% OPTIONS (name-value pairs)
%   'permutations'      : bootstrap replicates (default 5000)
%   'clusteralpha'      : two-sided cluster-forming threshold (default 0.05)
%   'randomseed'        : reproducible seed (default 1)
%   'time'              : optional time vector in seconds (default [])
%
% OUTPUT DIMENSIONS: time x selected networks x conditions x group pairs
%   significant         : samples belonging to significant clusters
%   corrected_p         : cluster-level p repeated over each observed cluster;
%                         NaN outside observed clusters (NOT pointwise p)
%   critical_p          : alpha (not an FDR critical p)
%   CLUSTER_STATS       : t statistics, clusters, and bootstrap null maxima
%                         Clusters{network,condition,pair} is a table of ONLY
%                         significant clusters, chronological and numbered 1,2,...
%                         NullMaximumClusterMass is replicates x networks x
%                         conditions x group pairs.
%
% CORRECTION FAMILY: all time samples separately for EACH network, condition,
% and group pair. There is NO correction across networks, conditions, or pairs.
%
% RESAMPLING: each group's participant time series is centered on its own
% observed group mean (imposing equal means of zero under the null). Each
% participant's entire centered residual time series is multiplied by one
% random +/-1 multiplier. Welch t statistics are recomputed in each draw.
% This is an APPROXIMATE heteroscedasticity-aware wild bootstrap, NOT an
% exact permutation test. It requires independent participants and is not
% guaranteed to control error in small/unbalanced groups; validate with
% null simulations before confirmatory use. In particular, groups with only
% 2 participants can have unstable bootstrap variances.
%
% A significant cluster supports cluster-level inference, not exact onset
% or individual-sample significance. Requires Statistics Toolbox (tinv).

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
if nTime < 2 || nNetworks < 1 || nConditions < 1 || nParticipants < 4
    error('Need >=2 time points, >=1 network, >=1 condition, >=4 participants.');
end
if ~isnumeric(group_labels) || numel(group_labels) ~= nParticipants || ...
        any(~isfinite(group_labels(:))) || any(group_labels(:) < 1) || ...
        any(group_labels(:) ~= floor(group_labels(:)))
    error('group_labels must contain one positive integer per participant.');
end
group_labels = group_labels(:);
if ~isnumeric(group_pairs) || isempty(group_pairs) || ...
        size(group_pairs,2) ~= 2 || any(~isfinite(group_pairs(:))) || ...
        any(group_pairs(:) < 1) || any(group_pairs(:) ~= floor(group_pairs(:))) || ...
        any(group_pairs(:,1) == group_pairs(:,2)) || ...
        any(~ismember(group_pairs(:),unique(group_labels))) || ...
        size(unique(sort(group_pairs,2),'rows'),1) ~= size(group_pairs,1)
    error('group_pairs must contain distinct valid nonduplicated group pairs.');
end
if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 0 || alpha >= 1 || ...
        ~isscalar(opts.clusteralpha) || ~isfinite(opts.clusteralpha) || ...
        opts.clusteralpha <= 0 || opts.clusteralpha >= 1
    error('alpha and clusteralpha must lie strictly between 0 and 1.');
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
        error('time must be a finite, strictly increasing vector with nTime values.');
    end
    opts.time = opts.time(:);
end

nPairs = size(group_pairs,1);
group_indices = cell(nPairs,2);
for pair = 1:nPairs
    for g = 1:2
        group_indices{pair,g} = find(group_labels == group_pairs(pair,g));
        if numel(group_indices{pair,g}) < 2
            error('Each group in a requested pair needs at least 2 participants.');
        end
    end
end

significant = false(nTime,nNetworks,nConditions,nPairs);
corrected_p = nan(nTime,nNetworks,nConditions,nPairs);
critical_p = alpha;
t_statistics = zeros(nTime,nNetworks,nConditions,nPairs);
cluster_labels = zeros(nTime,nNetworks,nConditions,nPairs);
cluster_tables = cell(nNetworks,nConditions,nPairs);
null_maximum_mass = zeros(opts.permutations,nNetworks,nConditions,nPairs);
cluster_threshold = nan(nTime,nNetworks,nConditions,nPairs);

previous_rng = rng;
cleanup_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(opts.randomseed,'twister');

for net = 1:nNetworks
    for condition = 1:nConditions
        for pair = 1:nPairs
            x = reshape(network_time_series(:,net,condition,group_indices{pair,1}), ...
                nTime,[]);
            y = reshape(network_time_series(:,net,condition,group_indices{pair,2}), ...
                nTime,[]);
            [observed_t,observed_df] = welch_t(x,y);
            observed_threshold = tinv(1-opts.clusteralpha/2,observed_df);
            t_statistics(:,net,condition,pair) = observed_t;
            cluster_threshold(:,net,condition,pair) = observed_threshold;
            [labels,masses] = temporal_clusters(observed_t,observed_threshold);
            cluster_labels(:,net,condition,pair) = labels;

            % Center residuals separately within each group to impose the null.
            residual_x = bsxfun(@minus,x,mean(x,2));
            residual_y = bsxfun(@minus,y,mean(y,2));
            for perm = 1:opts.permutations
                signs_x = 2*(rand(1,size(x,2)) > 0.5)-1;
                signs_y = 2*(rand(1,size(y,2)) > 0.5)-1;
                xb = bsxfun(@times,residual_x,signs_x);
                yb = bsxfun(@times,residual_y,signs_y);
                [boot_t,boot_df] = welch_t(xb,yb);
                boot_threshold = tinv(1-opts.clusteralpha/2,boot_df);
                [~,boot_masses] = temporal_clusters(boot_t,boot_threshold);
                if ~isempty(boot_masses)
                    null_maximum_mass(perm,net,condition,pair) = max(boot_masses);
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
                current_t = observed_t(members);
                mean_stat(c) = mean(current_t);
                [~,ix] = max(abs(current_t));
                peak_stat(c) = current_t(ix);
                cluster_p(c) = (1+sum(null_maximum_mass(:,net,condition,pair) ...
                    >= masses(c))) / (opts.permutations+1);
                corrected_p(members,net,condition,pair) = cluster_p(c);
                significant(members,net,condition,pair) = cluster_p(c) <= alpha;
            end
            keep = find(cluster_p <= alpha);
            if isempty(opts.time)
                [~,order] = sort(starts(keep));
            else
                [~,order] = sort(opts.time(starts(keep)));
            end
            keep = keep(order);
            cluster_number = (1:numel(keep))';
            if isempty(opts.time)
                cluster_tables{net,condition,pair} = table(cluster_number, ...
                    starts(keep),ends(keep),mean_stat(keep),peak_stat(keep), ...
                    masses(keep),cluster_p(keep), ...
                    'VariableNames',{'Cluster','StartSample','EndSample', ...
                    'MeanStatistic','PeakStatistic','Mass','CorrectedP'});
            else
                cluster_tables{net,condition,pair} = table(cluster_number, ...
                    opts.time(starts(keep)),opts.time(ends(keep)), ...
                    mean_stat(keep),peak_stat(keep),masses(keep),cluster_p(keep), ...
                    'VariableNames',{'Cluster','StartTime','EndTime', ...
                    'MeanStatistic','PeakStatistic','Mass','CorrectedP'});
            end
        end
    end
end

CLUSTER_STATS.Settings.Alpha = alpha;
CLUSTER_STATS.Settings.ClusterFormingAlpha = opts.clusteralpha;
CLUSTER_STATS.Settings.Permutations = opts.permutations;
CLUSTER_STATS.Settings.RandomSeed = opts.randomseed;
CLUSTER_STATS.Settings.Time = opts.time;
CLUSTER_STATS.Settings.GroupPairs = group_pairs;
CLUSTER_STATS.Settings.GroupLabels = group_labels;
CLUSTER_STATS.Settings.CorrectionFamily = ...
    'All time points, separately per network, condition and group pair';
CLUSTER_STATS.Settings.Resampling = ...
    'Within-group centered-residual Rademacher wild bootstrap';
CLUSTER_STATS.Settings.InferenceCaveat = ...
    'Approximate cluster-level inference; validate with null simulations, especially for small or unequal groups';
CLUSTER_STATS.TStatistics = t_statistics;
CLUSTER_STATS.ClusterFormingThreshold = cluster_threshold;
CLUSTER_STATS.ClusterLabels = cluster_labels;
CLUSTER_STATS.Clusters = cluster_tables;
CLUSTER_STATS.NullMaximumClusterMass = null_maximum_mass;
CLUSTER_STATS.ClusterCorrectedP = corrected_p;
CLUSTER_STATS.Significant = significant;
end

function [t,df] = welch_t(x,y)
    n1 = size(x,2);
    n2 = size(y,2);
    v1 = var(x,0,2)/n1;
    v2 = var(y,0,2)/n2;
    denominator = sqrt(v1+v2);
    t = (mean(x,2)-mean(y,2)) ./ denominator;
    df = (v1+v2).^2 ./ (v1.^2/(n1-1) + v2.^2/(n2-1));
    % Degenerate bootstrap samples cannot form valid suprathreshold clusters.
    t(~isfinite(t)) = 0;
    df(~isfinite(df) | df <= 0) = n1+n2-2;
end

function [labels,masses] = temporal_clusters(t,threshold)
    labels = zeros(numel(t),1);
    masses = zeros(0,1);
    for direction = [1 -1]
        active = (direction*t(:)) >= threshold(:);
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
