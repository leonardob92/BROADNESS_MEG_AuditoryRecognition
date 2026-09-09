% ========================================================================
%  BROADNESS STATISTICS EXAMPLE 3: PARTICIPANT GROUPS AND CONDITIONS
% ========================================================================
%
%  This script demonstrates comparisons between two or more independent
%  participant groups, including the condition x group interaction.
%
%  IMPORTANT DATA ORGANIZATION:
%  Concatenate every participant from every group along the fourth input
%  dimension and estimate BROADNESS only once. Keep a separate group label
%  for each participant. A common decomposition is necessary because PCA
%  or ICA estimated separately in each group would not guarantee matching
%  component order, sign, or scale.
%
%  Statistics and Machine Learning Toolbox is required.
% ========================================================================

%% 0) STARTUP AND INPUT CHECKS

toolbox_path = fileparts(mfilename('fullpath'));
addpath(toolbox_path)
BROADNESS_Startup(toolbox_path);

if ~exist('BROADNESS', 'var') || ...
        ~isfield(BROADNESS, 'TimeSeries_BrainNetworks') || ...
        ~isfield(BROADNESS, 'Time')
    error('Provide a participant-level BROADNESS structure in the workspace.');
end
required_functions = {'ttest2','fitlme'};
for function_index = 1:length(required_functions)
    if exist(required_functions{function_index}, 'file') ~= 2
        error(['This example requires Statistics and Machine Learning Toolbox. ' ...
            'Missing function: ' required_functions{function_index} '.']);
    end
end

network_time_series = BROADNESS.TimeSeries_BrainNetworks;
time = BROADNESS.Time(:);
[number_timepoints, number_networks, number_conditions, number_participants] = ...
    size(network_time_series);
if number_timepoints ~= length(time) || number_participants < 4
    error('Participant-level data and matching time information are required.');
end

%% 1) USER SETTINGS

% Define one group label per participant in exactly the same order as the
% fourth dimension of BROADNESS.TimeSeries_BrainNetworks. Numeric labels,
% strings, or categorical values are accepted.
if ~exist('group_labels', 'var')
    group_labels = []; %example: [ones(1,20) 2*ones(1,22)]
end
if ~exist('behaviour', 'var')
    behaviour = []; %optional: one numeric value per participant
end

selected_networks = 1:min(2,number_networks);
condition_labels = arrayfun(@(condition) ['Condition ' num2str(condition)], ...
    1:number_conditions, 'UniformOutput', false);
summary_time_window = [0.35 0.75]; %choose this interval a priori
alpha = 0.05;
minimum_significant_samples = 1;

figure_mode = 'show';
figure_layout = 'individual';
output_path = fullfile(fileparts(toolbox_path), 'Output');

if isempty(group_labels) || numel(group_labels) ~= number_participants
    error(['Provide one entry in "group_labels" for every participant in ' ...
        'the fourth time-series dimension.']);
end
if length(condition_labels) ~= number_conditions
    error('Provide one entry in "condition_labels" for every condition.');
end
summary_indices = time >= summary_time_window(1) & time <= summary_time_window(2);
if ~any(summary_indices)
    error('"summary_time_window" does not overlap the available time vector.');
end

categorical_groups = categorical(group_labels(:));
group_categories = categories(categorical_groups);
number_groups = length(group_categories);
if number_groups < 2
    error('At least two participant groups are required.');
end
group_indices = cell(number_groups,1);
group_names = cell(number_groups,1);
for group = 1:number_groups
    group_indices{group} = find(categorical_groups == group_categories{group});
    group_names{group} = char(group_categories{group});
    if length(group_indices{group}) < 2
        error('Every group must contain at least two participants.');
    end
end

%% 2) INDEPENDENT GROUP TESTS FOR EACH CONDITION

% Welch t-tests are used because they do not assume equal group variances.
% ranksum can replace them time-point by time-point when a nonparametric
% group contrast is required; the same FDR family should then be retained.

group_pairs = nchoosek(1:number_groups,2);
number_group_pairs = size(group_pairs,1);
p_values = nan(number_timepoints,length(selected_networks), ...
    number_conditions,number_group_pairs);
t_statistics = nan(size(p_values));
hedges_g = nan(size(p_values));

for network_index = 1:length(selected_networks)
    network = selected_networks(network_index);
    for condition = 1:number_conditions
        for pair = 1:number_group_pairs
            group_1 = group_pairs(pair,1);
            group_2 = group_pairs(pair,2);
            values_1 = reshape(network_time_series(:,network,condition, ...
                group_indices{group_1}),number_timepoints,[]);
            values_2 = reshape(network_time_series(:,network,condition, ...
                group_indices{group_2}),number_timepoints,[]);
            [~,p_values(:,network_index,condition,pair),~,test_statistics] = ...
                ttest2(values_1,values_2,'Dim',2,'Vartype','unequal');
            t_statistics(:,network_index,condition,pair) = test_statistics.tstat;

            number_1 = sum(~isnan(values_1),2);
            number_2 = sum(~isnan(values_2),2);
            variance_1 = var(values_1,0,2,'omitnan');
            variance_2 = var(values_2,0,2,'omitnan');
            pooled_sd = sqrt(((number_1-1).*variance_1 + (number_2-1).*variance_2) ./ ...
                (number_1+number_2-2));
            cohen_d = (mean(values_1,2,'omitnan') - mean(values_2,2,'omitnan')) ./ pooled_sd;
            small_sample_correction = 1 - 3 ./ (4*(number_1+number_2)-9);
            hedges_g(:,network_index,condition,pair) = ...
                small_sample_correction .* cohen_d;
        end
    end
end
hedges_g(~isfinite(hedges_g)) = NaN;

[significant,adjusted_p,critical_p] = BROADNESS_FDRCorrection(p_values,alpha);

group_contrast_labels = cell(number_conditions,number_group_pairs);
for condition = 1:number_conditions
    for pair = 1:number_group_pairs
        group_contrast_labels{condition,pair} = ...
            [condition_labels{condition} ': ' ...
            group_names{group_pairs(pair,1)} ' vs ' ...
            group_names{group_pairs(pair,2)}];
    end
end

%% 3) CONDITION x GROUP MIXED MODEL IN A PREDEFINED TIME WINDOW

window_models = cell(length(selected_networks),1);
window_tables = cell(length(selected_networks),1);
for network_index = 1:length(selected_networks)
    network = selected_networks(network_index);
    condition_means = reshape(mean(network_time_series(summary_indices,network,:,:), ...
        1,'omitnan'),number_conditions,number_participants)';
    response = reshape(condition_means',[],1);
    participant = categorical(repelem((1:number_participants)',number_conditions));
    condition = categorical(repmat((1:number_conditions)',number_participants,1));
    group = repelem(categorical_groups,number_conditions);
    window_tables{network_index} = table(response,participant,condition,group, ...
        'VariableNames', {'Response','Participant','Condition','Group'});
    window_models{network_index} = fitlme(window_tables{network_index}, ...
        'Response ~ Condition*Group + (1|Participant)');
end

%% 4) OPTIONAL BEHAVIOURAL MIXED-EFFECTS MODEL

behaviour_models = cell(length(selected_networks),1);
if ~isempty(behaviour)
    if ~isnumeric(behaviour) || numel(behaviour) ~= number_participants
        error('"behaviour" must contain one numeric value per participant.');
    end
    for network_index = 1:length(selected_networks)
        behaviour_long = repelem(behaviour(:),number_conditions);
        behaviour_table = window_tables{network_index};
        behaviour_table.Behaviour = behaviour_long;
        behaviour_models{network_index} = fitlme(behaviour_table, ...
            'Response ~ Condition*Group + Behaviour + (1|Participant)');
    end
end

%% 5) PLOT THE GROUP CONTRASTS

number_contrasts = number_conditions * number_group_pairs;
significance_colors = lines(max(1,number_contrasts));
condition_group_colors = cell(1,number_conditions);
all_colors = lines(number_conditions*number_groups);
for condition = 1:number_conditions
    color_indices = condition:number_conditions:(number_conditions*number_groups);
    condition_group_colors{condition} = all_colors(color_indices,:);
end

statistical_figures = cell(length(selected_networks),1);
for network_index = 1:length(selected_networks)
    significant_windows = cell(number_contrasts,1);
    contrast = 0;
    for condition = 1:number_conditions
        for pair = 1:number_group_pairs
            contrast = contrast + 1;
            significant_windows{contrast} = BROADNESS_SignificantSamplesToWindows( ...
                significant(:,network_index,condition,pair),time,minimum_significant_samples);
        end
    end
    statistical_figures{network_index} = BROADNESS_Plot_ActivationTimeseries( ...
        network_time_series(:,selected_networks(network_index),:,:),time, ...
        'Groups',group_indices,'GroupLabels',group_names, ...
        'ConditionLabels',condition_labels, ...
        'ConditionXGroupColors',condition_group_colors, ...
        'SignificantWindows',significant_windows, ...
        'SignificanceColors',significance_colors, ...
        'SignificanceStyle','line', ...
        'FigureMode',figure_mode,'FigureLayout',figure_layout, ...
        'OutputPath',output_path,'FigurePrefix', ...
        ['Groups_BN' num2str(selected_networks(network_index))]);
end

%% 6) COLLECT THE STATISTICAL OUTPUTS

GROUP_STATS.SelectedNetworks = selected_networks;
GROUP_STATS.GroupPairs = group_pairs;
GROUP_STATS.GroupNames = group_names;
GROUP_STATS.ContrastLabels = group_contrast_labels;
GROUP_STATS.Time = time;
GROUP_STATS.PValues = p_values;
GROUP_STATS.AdjustedPValues = adjusted_p;
GROUP_STATS.Significant = significant;
GROUP_STATS.FDRCriticalP = critical_p;
GROUP_STATS.TStatistics = t_statistics;
GROUP_STATS.HedgesG = hedges_g;
GROUP_STATS.SummaryTimeWindow = summary_time_window;
GROUP_STATS.WindowTables = window_tables;
GROUP_STATS.WindowModels = window_models;
GROUP_STATS.BehaviourModels = behaviour_models;
GROUP_STATS.Figures = statistical_figures;

disp('Group statistical example completed. Results are stored in GROUP_STATS.');
