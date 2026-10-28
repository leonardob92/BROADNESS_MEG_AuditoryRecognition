% ========================================================================
%  BROADNESS STATISTICS EXAMPLE 2: INDEPENDENT PARTICIPANT GROUPS
% ========================================================================
%
%  Compare two or more independent participant groups within each condition.
%  All participants are concatenated along the fourth dimension BEFORE
%  running one common BROADNESS network estimation.
%  Network time series: time x networks x conditions x participants.
%  The participant (not the time sample) is the independent unit.
%  Optional window-level models test condition x group and behaviour.
%
%  Requires Statistics and Machine Learning Toolbox (ttest2, fitlme).
% ========================================================================

%% 0) STARTUP

% Simply download the BROADNESS Toolbox folder and place it in your working directory,
% making sure not to alter the structure of its functions, subfolders, or files.

clear 
close all
clc

% Setup directories relative to this example script
project_path = '/Users/au550322/Documents/GitHub/BROADNESS_MEG_AuditoryRecognition/BROADNESS_Toolbox';
data_path = '/Users/au550322/Documents/GitHub/BROADNESS_MEG_AuditoryRecognition/Data';
output_path = [project_path '/Output/Stats'];
addpath(project_path)
BROADNESS_Startup(project_path);

%% 1) PERFORM BROADNESS (ONLY ESSENTIAL INPUTS)

%%% ------------------- USER SETTINGS ------------------- %%%

% loading a few participants and concatenating them
list = dir(fullfile(data_path, 'SUBJ*.mat'));
data = [];
for subi = 1:length(list) %over participants
    load(fullfile(list(subi).folder, list(subi).name)) %loading data for each participant
    data = cat(4,data,Data); %concatenating data
    disp(subi)
end
data(:,777:end,:,:) = []; %this is simply because we have too many data points in the single-participant data in this example
%loading time
load(fullfile(data_path, 'DataReduced_AveragedOverParticipants_Example.mat'),'time');

%%% ------------------ COMPUTATION --------------------- %%%

% Run BROADNESS network estimation (default parameters)
BROADNESS = BROADNESS_NetworkEstimation(data, time);

% Extracting relevant information from BROADNESS_NetworkEstimation
network_time_series = BROADNESS.TimeSeries_BrainNetworks;
time = BROADNESS.Time(:);
[number_timepoints, number_networks, number_conditions, number_participants] = ...
    size(network_time_series);
if number_timepoints ~= length(time) || number_participants < 2 || number_conditions < 2
    error(['This example requires matching time information, at least two ' ...
        'conditions, and participant-level data from at least two participants.']);
end


%% 2) USER SETTINGS

selected_networks = [1 2];
condition_labels = {'Memorized','NewT1','NewT2','NewT3','NewT4'};
alpha = 0.05;
minimum_significant_samples = 1;

% IMPORTANT: One group label for EACH participant, in the order of "list".
% Example for 10 participants: first 5 in Control, next 5 in Patient.
% For 3 groups, use e.g. [ones(1,5) 2*ones(1,5) 3*ones(1,5)]
% and provide three group names (matching the participant file order).
% REPLACE this example with your actual participant assignments.
group_labels = [ones(1,20) 2*ones(1,15) 3*ones(1,8)];
group_names = {'Control','Patient','s'}; %group 1 and group 2

% Predefine the time window for the condition x group interaction.
summary_time_window = [0.35 0.75]; %seconds
run_interaction = 'on'; %'on' or 'off'

% Optional: one independent behavioural value per participant, in list order.
% Leave empty to skip. Do not derive behaviour from the tested MEG data.
behaviour = [];

figure_mode = 'show'; %'off', 'show', 'save', or 'both'
figure_layout = 'individual';

% Check the settings before starting the statistics.
if numel(condition_labels) ~= number_conditions
    error('Provide one condition label for each condition.');
end
if any(selected_networks < 1) || any(selected_networks > number_networks) || ...
        any(mod(selected_networks,1) ~= 0)
    error('selected_networks contains an invalid network index.');
end
number_groups = numel(group_names);
if number_groups < 2 || numel(group_labels) ~= number_participants || ...
        any(~ismember(group_labels,1:number_groups))
    error('Provide one group label (1, 2, ...) per participant and matching group names.');
end
for group = 1:number_groups
    if sum(group_labels == group) < 2
        error('Every group must contain at least two participants.');
    end
end
if ~isempty(behaviour) && (~isnumeric(behaviour) || ...
        numel(behaviour) ~= number_participants || any(~isfinite(behaviour(:))))
    error('behaviour must contain one finite numeric value per participant.');
end
if ~ismember(lower(run_interaction),{'on','off'})
    error('run_interaction must be ''on'' or ''off''.');
end
summary_indices = time >= summary_time_window(1) & ...
    time <= summary_time_window(2);
if ~any(summary_indices)
    error('summary_time_window does not overlap the time vector.');
end

group_indices = cell(number_groups,1);
for group = 1:number_groups
    group_indices{group} = find(group_labels == group);
end

%% 3) INDEPENDENT GROUP TESTS FOR EACH CONDITION

% At every time point, condition and network, compare independent groups.
% With 2 groups: Welch independent t-test.
% With 3+ groups: Welch one-way ANOVA (overall group effect), plus Welch
% t-tests for every group pair. Welch tests allow unequal group variances.
% Pairwise effect size: Hedges' g (first group minus second group).
% Overall and pairwise p-values use SEPARATE FDR families, each including
% all selected time points, networks and conditions (and all group pairs
% for pairwise tests). Pairwise tests are not gated on ANOVA significance.

group_pairs = nchoosek(1:number_groups,2);
number_group_pairs = size(group_pairs,1);
number_selected_networks = length(selected_networks);

% Dimensions: time x selected networks x conditions x group pairs.
p_values = nan(number_timepoints,number_selected_networks, ...
    number_conditions,number_group_pairs);
t_statistics = nan(size(p_values));
hedges_g = nan(size(p_values));

% Overall ANOVA is only needed when there are 3 or more groups.
overall_p = nan(number_timepoints,number_selected_networks,number_conditions);
overall_F = nan(size(overall_p));

for network_index = 1:number_selected_networks
    network = selected_networks(network_index);
    for condition = 1:number_conditions
        group_values = cell(number_groups,1);
        for group = 1:number_groups
            group_values{group} = reshape( ...
                network_time_series(:,network,condition,group_indices{group}), ...
                number_timepoints,[]);
        end

        % Welch one-way ANOVA, calculated separately at each time point.
        % This is the unequal-variance alternative to ordinary ANOVA.
        if number_groups > 2
            for sample = 1:number_timepoints
                n = zeros(number_groups,1);
                m = zeros(number_groups,1);
                v = zeros(number_groups,1);
                for group = 1:number_groups
                    x = group_values{group}(sample,:);
                    x = x(~isnan(x));
                    n(group) = numel(x);
                    if n(group) >= 2
                        m(group) = mean(x);
                        v(group) = var(x,0);
                    end
                end
                if any(n < 2) || any(v <= 0) || any(~isfinite(v))
                    continue
                end
                w = n ./ v;
                wsum = sum(w);
                weighted_mean = sum(w.*m) / wsum;
                correction_sum = sum((1-w/wsum).^2 ./ (n-1));
                df1 = number_groups-1;
                df2 = (number_groups^2-1) / (3*correction_sum);
                F = (sum(w.*(m-weighted_mean).^2)/df1) / ...
                    (1 + 2*(number_groups-2)*correction_sum / ...
                    (number_groups^2-1));
                overall_F(sample,network_index,condition) = F;
                overall_p(sample,network_index,condition) = fcdf(F,df1,df2,'upper');
            end
        end

        % Pairwise Welch t-tests and Hedges' g.
        for pair = 1:number_group_pairs
            values_1 = group_values{group_pairs(pair,1)};
            values_2 = group_values{group_pairs(pair,2)};
            [~,p_values(:,network_index,condition,pair),~,test_statistics] = ...
                ttest2(values_1,values_2,'Dim',2,'Vartype','unequal');
            t_statistics(:,network_index,condition,pair) = test_statistics.tstat;

            n1 = sum(~isnan(values_1),2);
            n2 = sum(~isnan(values_2),2);
            variance_1 = var(values_1,0,2,'omitnan');
            variance_2 = var(values_2,0,2,'omitnan');
            pooled_sd = sqrt(((n1-1).*variance_1 + (n2-1).*variance_2) ./ ...
                (n1+n2-2));
            cohen_d = (mean(values_1,2,'omitnan') - ...
                mean(values_2,2,'omitnan')) ./ pooled_sd;
            correction = 1 - 3 ./ (4*(n1+n2)-9);
            hedges_g(:,network_index,condition,pair) = correction .* cohen_d;
        end
    end
end
hedges_g(~isfinite(hedges_g)) = NaN;
[significant,adjusted_p,critical_p] = BROADNESS_FDRCorrection(p_values,alpha);

if number_groups > 2
    [overall_significant,overall_adjusted_p,overall_critical_p] = ...
        BROADNESS_FDRCorrection(overall_p,alpha);
else
    overall_significant = false(size(overall_p));
    overall_adjusted_p = nan(size(overall_p));
    overall_critical_p = NaN;
end

GROUP_STATS.Settings.SelectedNetworks = selected_networks;
GROUP_STATS.Settings.ConditionLabels = condition_labels;
GROUP_STATS.Settings.GroupNames = group_names;
GROUP_STATS.Settings.GroupLabels = group_labels;
GROUP_STATS.Settings.NumberGroups = number_groups;
GROUP_STATS.Settings.Time = time;
GROUP_STATS.Settings.Alpha = alpha;
GROUP_STATS.Settings.MinimumSignificantSamples = minimum_significant_samples;
GROUP_STATS.Settings.FigureMode = figure_mode;
GROUP_STATS.Settings.FigureLayout = figure_layout;
GROUP_STATS.BetweenGroupComparisons.GroupPairs = group_pairs;
GROUP_STATS.BetweenGroupComparisons.PValues = p_values;
GROUP_STATS.BetweenGroupComparisons.AdjustedPValues = adjusted_p;
GROUP_STATS.BetweenGroupComparisons.Significant = significant;
GROUP_STATS.BetweenGroupComparisons.FDRCriticalP = critical_p;
GROUP_STATS.BetweenGroupComparisons.FDRMethod = 'Benjamini-Hochberg';
GROUP_STATS.BetweenGroupComparisons.FDRFamily = ...
    'All selected time-points, networks, conditions and group pairs';
GROUP_STATS.BetweenGroupComparisons.TStatistics = t_statistics;
GROUP_STATS.BetweenGroupComparisons.HedgesG = hedges_g;
GROUP_STATS.BetweenGroupComparisons.OverallANOVA.Enabled = number_groups > 2;
GROUP_STATS.BetweenGroupComparisons.OverallANOVA.FStatistics = overall_F;
GROUP_STATS.BetweenGroupComparisons.OverallANOVA.PValues = overall_p;
GROUP_STATS.BetweenGroupComparisons.OverallANOVA.AdjustedPValues = overall_adjusted_p;
GROUP_STATS.BetweenGroupComparisons.OverallANOVA.Significant = overall_significant;
GROUP_STATS.BetweenGroupComparisons.OverallANOVA.FDRCriticalP = overall_critical_p;
GROUP_STATS.BetweenGroupComparisons.OverallANOVA.FDRMethod = 'Benjamini-Hochberg';
GROUP_STATS.BetweenGroupComparisons.OverallANOVA.FDRFamily = ...
    'All selected time-points, networks and conditions';

%% 4) PLOT THE FDR-SIGNIFICANT WINDOWS

% With 2 groups: one numbered significance line per condition.
% With 3+ groups: A1, A2, ... show the overall ANOVA for each condition;
% numbered lines show pairwise group comparisons within each condition.
% These are FDR-significant samples, not cluster-corrected intervals.

number_pairwise_contrasts = number_conditions*number_group_pairs;
pairwise_codes = arrayfun(@num2str,1:number_pairwise_contrasts, ...
    'UniformOutput',false)';
pairwise_descriptions = cell(number_pairwise_contrasts,1);
contrast = 0;
for condition = 1:number_conditions
    for pair = 1:number_group_pairs
        contrast = contrast + 1;
        pairwise_descriptions{contrast} = [condition_labels{condition} ': ' ...
            group_names{group_pairs(pair,1)} ' vs ' ...
            group_names{group_pairs(pair,2)}];
    end
end

if number_groups > 2
    anova_codes = arrayfun(@(c) ['A' num2str(c)],1:number_conditions, ...
        'UniformOutput',false)';
    anova_descriptions = cell(number_conditions,1);
    for condition = 1:number_conditions
        anova_descriptions{condition} = ...
            [condition_labels{condition} ': overall group ANOVA'];
    end
    significance_labels = [anova_codes;pairwise_codes];
    significance_descriptions = [anova_descriptions;pairwise_descriptions];
else
    significance_labels = pairwise_codes;
    significance_descriptions = pairwise_descriptions;
end
significance_colors = lines(numel(significance_labels));
all_colors = lines(number_conditions*number_groups);
condition_group_colors = cell(1,number_conditions);
for condition = 1:number_conditions
    condition_group_colors{condition} = ...
        all_colors(condition:number_conditions:end,:);
end
significance_key = table(significance_labels,significance_descriptions, ...
    'VariableNames',{'Line','StatisticalComparison'});
disp('Significance-line key:')
disp(significance_key)

statistical_figures = cell(number_selected_networks,1);
for network_index = 1:number_selected_networks
    significant_windows = cell(numel(significance_labels),1);
    contrast = 0;
    if number_groups > 2
        for condition = 1:number_conditions
            contrast = contrast + 1;
            significant_windows{contrast} = BROADNESS_SignificantSamplesToWindows( ...
                overall_significant(:,network_index,condition),time, ...
                minimum_significant_samples);
        end
    end
    for condition = 1:number_conditions
        for pair = 1:number_group_pairs
            contrast = contrast + 1;
            significant_windows{contrast} = BROADNESS_SignificantSamplesToWindows( ...
                significant(:,network_index,condition,pair),time, ...
                minimum_significant_samples);
        end
    end
    statistical_figures{network_index} = BROADNESS_Plot_ActivationTimeseries( ...
        network_time_series(:,selected_networks(network_index),:,:),time, ...
        'Groups',group_indices,'GroupLabels',group_names, ...
        'ConditionLabels',condition_labels, ...
        'ConditionXGroupColors',condition_group_colors, ...
        'SignificantWindows',significant_windows, ...
        'SignificanceColors',significance_colors, ...
        'SignificanceLabels',significance_labels, ...
        'SignificanceStyle','line', ...
        'SignificanceLineOffset',0.025, ...
        'FigureMode',figure_mode,'FigureLayout',figure_layout, ...
        'OutputPath',output_path,'FigurePrefix', ...
        ['BetweenGroups_BN' num2str(selected_networks(network_index))]);
end
GROUP_STATS.Figures.SignificanceKey = significance_key;
GROUP_STATS.Figures.GroupComparisons = statistical_figures;

%% 5) CONDITION x GROUP INTERACTION IN A PREDEFINED TIME WINDOW

% This asks a DIFFERENT question from Section 3:
% "Does the effect of condition depend on which group a participant belongs to?"
% Average each participant's network response over a PREDEFINED time window.
% Then fit a mixed-effects model: Condition and Group are fixed effects;
% Participant has a random intercept because conditions repeat within person.
% This is ONE window-level model per network, not a time-resolved test.
% Inspect the Condition:Group term(s) with anova(model).

window_models = cell(length(selected_networks),1);
window_anova = cell(length(selected_networks),1);
window_tables = cell(length(selected_networks),1);
if strcmpi(run_interaction,'on')
    for network_index = 1:length(selected_networks)
        network = selected_networks(network_index);
        condition_means = reshape(mean( ...
            network_time_series(summary_indices,network,:,:),1,'omitnan'), ...
            number_conditions,number_participants)';
        response = reshape(condition_means',[],1);
        participant = categorical(repelem((1:number_participants)',number_conditions));
        condition = categorical(repmat((1:number_conditions)',number_participants,1));
        group = categorical(repelem(group_labels(:),number_conditions));
        window_tables{network_index} = table(response,participant,condition,group, ...
            'VariableNames',{'Response','Participant','Condition','Group'});
        window_models{network_index} = fitlme(window_tables{network_index}, ...
            'Response ~ Condition*Group + (1|Participant)');
        window_anova{network_index} = anova(window_models{network_index});
        disp(['Condition x group model: Network ' num2str(network)])
        disp(window_anova{network_index})
    end
end
GROUP_STATS.ConditionGroupMixedModel.Enabled = strcmpi(run_interaction,'on');
GROUP_STATS.ConditionGroupMixedModel.TimeWindow = summary_time_window;
GROUP_STATS.ConditionGroupMixedModel.Tables = window_tables;
GROUP_STATS.ConditionGroupMixedModel.Models = window_models;
GROUP_STATS.ConditionGroupMixedModel.ANOVA = window_anova;

%% 6) OPTIONAL BEHAVIOURAL ASSOCIATION

% If behaviour is provided, fit a window-level mixed-effects model that
% includes behaviour as an additional participant-level predictor.
% This estimates an association adjusted for condition and group; it does
% NOT test whether the behaviour association differs between groups.
% It does not establish causality. Leave behaviour = [] to skip.

% HOW TO INTERPRET THE RESULTS:
% Open GROUP_STATS.Behaviour.Models{1}.Coefficients for Network 1,
% or Models{2}.Coefficients for Network 2.
%
% The Coefficients table contains 8 columns:
%   Name     : Model coefficient (Condition, Group, Behaviour, etc.).
%   Estimate : Estimated effect size and direction.
%   SE       : Standard error of the estimate.
%   tStat    : t-statistic testing whether the coefficient differs from zero.
%   DF       : Degrees of freedom for the t-test.
%   pValue   : Statistical significance (p < 0.05, uncorrected).
%   Lower    : Lower bound of the 95% confidence interval.
%   Upper    : Upper bound of the 95% confidence interval.
%
% For behavioural associations, find the row named 'Behaviour':
%   Estimate > 0 : Higher behaviour is associated with higher activation.
%   Estimate < 0 : Higher behaviour is associated with lower activation.
%   pValue < 0.05: Evidence of an association (uncorrected).
%
% The model accounts for Condition, Group and repeated observations
% from the same participant.

behaviour_models = cell(length(selected_networks),1);
if ~isempty(behaviour)
    if ~strcmpi(run_interaction,'on')
        error('Turn run_interaction on to prepare the window-level tables.');
    end
    for network_index = 1:length(selected_networks)
        behaviour_table = window_tables{network_index};
        behaviour_table.Behaviour = repelem(behaviour(:),number_conditions);
        behaviour_models{network_index} = fitlme(behaviour_table, ...
            'Response ~ Condition*Group + Behaviour + (1|Participant)');
    end
end
GROUP_STATS.Behaviour.Enabled = ~isempty(behaviour);
GROUP_STATS.Behaviour.Values = behaviour;
GROUP_STATS.Behaviour.Models = behaviour_models;

disp('Between-group statistical example completed. Results are stored in GROUP_STATS.');

%% 7) GROUP DIFFERENCES IN SPATIAL ACTIVATION PATTERNS

% Compare participant-specific spatial maps from the SAME BROADNESS
% decomposition, separately for each experimental condition.
% With 2 groups: independent Welch t-tests; with 3+ groups: additionally
% an overall Welch ANOVA. Spatial cluster correction uses an approximate
% heteroscedasticity-aware wild bootstrap (not paired sign flips).
% The function also compares six participant-level spatial summaries.
% Cluster significance applies to whole clusters, not individual sources.
% Start with 100 permutations to check the workflow; for final inference
% use at least 5000 and validate assumptions and bootstrap calibration.

selected_spatial_networks = 1;
spatial_permutations = 100; %exploration only; use >=5000 for final inference
spatial_map_normalization = 'none'; %'none' or 'rms'

load(fullfile(project_path,'BROADNESS_External', ...
    'MNI152_8mm_coord_dyi.mat'),'MNI8');

[SPATIAL_GROUP_STATS,SPATIAL_GROUP_FIGURES] = ...
    BROADNESS_SpatialPatternGroupStatistics(BROADNESS,group_labels, ...
    'principalcomps',selected_spatial_networks, ...
    'groupnames',group_names,'conditionnames',condition_labels, ...
    'mni_coords',MNI8,'mapnormalization',spatial_map_normalization, ...
    'permutations',spatial_permutations, ...
    'figuremode',figure_mode,'figurelayout','both', ...
    'OutputPath',output_path,'figureprefix','BetweenGroups', ...
    'clusteralpha',0.05,'alpha',alpha);

% Inspect SPATIAL_GROUP_STATS.Voxelwise.Pairwise.Clusters for corrected
% cluster results; SPATIAL_GROUP_STATS.Voxelwise.Omnibus for 3+ groups;
% SPATIAL_GROUP_STATS.SummaryTests for spatial summary comparisons.

%%


