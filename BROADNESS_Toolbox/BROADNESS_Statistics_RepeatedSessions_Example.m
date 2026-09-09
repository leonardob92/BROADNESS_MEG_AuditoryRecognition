% ========================================================================
%  BROADNESS STATISTICS EXAMPLE 2: CONDITIONS ACROSS REPEATED SESSIONS
% ========================================================================
%
%  This script demonstrates comparisons when the same participants complete
%  the same experimental conditions in two or more sessions.
%
%  IMPORTANT DATA ORGANIZATION:
%  Estimate the networks only once after concatenating condition x session
%  combinations along the third input dimension. For example:
%
%      data_combined = cat(3,session1_data,session2_data);
%      BROADNESS = BROADNESS_NetworkEstimation(data_combined,time);
%
%  This preserves one common component basis across sessions. Estimating
%  PCA or ICA independently in each session would not guarantee matching
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
required_functions = {'ttest','fitrm','ranova','fitlme'};
for function_index = 1:length(required_functions)
    if exist(required_functions{function_index}, 'file') ~= 2
        error(['This example requires Statistics and Machine Learning Toolbox. ' ...
            'Missing function: ' required_functions{function_index} '.']);
    end
end

combined_time_series = BROADNESS.TimeSeries_BrainNetworks;
time = BROADNESS.Time(:);
[number_timepoints, number_networks, number_combined_conditions, number_participants] = ...
    size(combined_time_series);
if number_timepoints ~= length(time) || number_participants < 2
    error('Participant-level data and matching time information are required.');
end

%% 1) USER SETTINGS

number_conditions = 2;
number_sessions = 2;
condition_labels = arrayfun(@(condition) ['Condition ' num2str(condition)], ...
    1:number_conditions, 'UniformOutput', false);
session_labels = arrayfun(@(session) ['Session ' num2str(session)], ...
    1:number_sessions, 'UniformOutput', false);
selected_networks = 1:min(2,number_networks);
summary_time_window = [0.35 0.75]; %choose this interval a priori
alpha = 0.05;
minimum_significant_samples = 1;
if ~exist('behaviour', 'var')
    behaviour = []; %optional: one participant-level value
end

figure_mode = 'show';
figure_layout = 'individual';
output_path = fullfile(fileparts(toolbox_path), 'Output');

if number_conditions * number_sessions ~= number_combined_conditions
    error(['The third time-series dimension must contain number_conditions ' ...
        'x number_sessions entries.']);
end
if length(condition_labels) ~= number_conditions || ...
        length(session_labels) ~= number_sessions
    error('Provide one label for every condition and session.');
end
if number_sessions < 2
    error('At least two repeated sessions are required.');
end
summary_indices = time >= summary_time_window(1) & time <= summary_time_window(2);
if ~any(summary_indices)
    error('"summary_time_window" does not overlap the available time vector.');
end

% cat(3,session1_data,session2_data,...) stores every condition from the
% first session, followed by every condition from the second session, etc.
session_time_series = reshape(combined_time_series,number_timepoints, ...
    number_networks,number_conditions,number_sessions,number_participants);

%% 2) PAIRED SESSION TESTS FOR EACH CONDITION

% Paired t-tests are used here. A time-point-wise signrank analysis can be
% substituted when a nonparametric planned contrast is required.

session_pairs = nchoosek(1:number_sessions,2);
number_session_pairs = size(session_pairs,1);
p_values = nan(number_timepoints,length(selected_networks), ...
    number_conditions,number_session_pairs);
t_statistics = nan(size(p_values));
cohen_dz = nan(size(p_values));

for network_index = 1:length(selected_networks)
    network = selected_networks(network_index);
    for condition = 1:number_conditions
        for pair = 1:number_session_pairs
            session_1 = session_pairs(pair,1);
            session_2 = session_pairs(pair,2);
            values_1 = reshape(session_time_series(:,network,condition,session_1,:), ...
                number_timepoints,number_participants);
            values_2 = reshape(session_time_series(:,network,condition,session_2,:), ...
                number_timepoints,number_participants);
            [~,p_values(:,network_index,condition,pair),~,test_statistics] = ...
                ttest(values_1,values_2,'Dim',2);
            t_statistics(:,network_index,condition,pair) = test_statistics.tstat;
            paired_differences = values_1 - values_2;
            cohen_dz(:,network_index,condition,pair) = ...
                mean(paired_differences,2,'omitnan') ./ ...
                std(paired_differences,0,2,'omitnan');
        end
    end
end
cohen_dz(~isfinite(cohen_dz)) = NaN;

[significant,adjusted_p,critical_p] = BROADNESS_FDRCorrection(p_values,alpha);

session_contrast_labels = cell(number_conditions,number_session_pairs);
for condition = 1:number_conditions
    for pair = 1:number_session_pairs
        session_contrast_labels{condition,pair} = ...
            [condition_labels{condition} ': ' ...
            session_labels{session_pairs(pair,1)} ' vs ' ...
            session_labels{session_pairs(pair,2)}];
    end
end

%% 3) CONDITION x SESSION REPEATED-MEASURES MODEL IN A TIME WINDOW

number_repeated_measures = number_conditions * number_sessions;
response_names = arrayfun(@(measure) ['Measure' num2str(measure)], ...
    1:number_repeated_measures, 'UniformOutput', false);
within_condition = categorical(repmat((1:number_conditions)',number_sessions,1));
within_session = categorical(kron((1:number_sessions)',ones(number_conditions,1)));
within_design = table(within_condition,within_session, ...
    'VariableNames', {'Condition','Session'});
window_models = cell(length(selected_networks),1);
window_anova = cell(length(selected_networks),1);
window_responses = cell(length(selected_networks),1);

for network_index = 1:length(selected_networks)
    network = selected_networks(network_index);
    response_matrix = reshape(mean(session_time_series(summary_indices,network,:,:,:), ...
        1,'omitnan'),number_repeated_measures,number_participants)';
    window_responses{network_index} = response_matrix;
    statistics_table = array2table(response_matrix,'VariableNames',response_names);
    model_formula = [response_names{1} '-' response_names{end} ' ~ 1'];
    window_models{network_index} = fitrm(statistics_table,model_formula, ...
        'WithinDesign',within_design);
    window_anova{network_index} = ranova(window_models{network_index}, ...
        'WithinModel','Condition*Session');
end

%% 4) OPTIONAL BEHAVIOURAL MIXED-EFFECTS MODEL

behaviour_models = cell(length(selected_networks),1);
if ~isempty(behaviour)
    if ~isnumeric(behaviour) || numel(behaviour) ~= number_participants
        error('"behaviour" must contain one numeric value per participant.');
    end
    for network_index = 1:length(selected_networks)
        response_matrix = window_responses{network_index};
        response = reshape(response_matrix',[],1);
        participant = categorical(repelem((1:number_participants)', ...
            number_repeated_measures));
        condition = repmat(within_condition,number_participants,1);
        session = repmat(within_session,number_participants,1);
        behaviour_long = repelem(behaviour(:),number_repeated_measures);
        behaviour_table = table(response,participant,condition,session,behaviour_long, ...
            'VariableNames', {'Response','Participant','Condition','Session','Behaviour'});
        behaviour_models{network_index} = fitlme(behaviour_table, ...
            'Response ~ Condition*Session + Behaviour + (1|Participant)');
    end
end

%% 5) PLOT THE SESSION CONTRASTS

combined_labels = cell(1,number_repeated_measures);
for session = 1:number_sessions
    for condition = 1:number_conditions
        combined_index = condition + (session-1) * number_conditions;
        combined_labels{combined_index} = ...
            [condition_labels{condition} ' — ' session_labels{session}];
    end
end

number_contrasts = number_conditions * number_session_pairs;
significance_colors = lines(max(1,number_contrasts));
statistical_figures = cell(length(selected_networks),1);
for network_index = 1:length(selected_networks)
    significant_windows = cell(number_contrasts,1);
    contrast = 0;
    for condition = 1:number_conditions
        for pair = 1:number_session_pairs
            contrast = contrast + 1;
            significant_windows{contrast} = BROADNESS_SignificantSamplesToWindows( ...
                significant(:,network_index,condition,pair),time,minimum_significant_samples);
        end
    end
    statistical_figures{network_index} = BROADNESS_Plot_ActivationTimeseries( ...
        combined_time_series(:,selected_networks(network_index),:,:),time, ...
        'ConditionLabels',combined_labels, ...
        'SignificantWindows',significant_windows, ...
        'SignificanceColors',significance_colors, ...
        'SignificanceStyle','line', ...
        'FigureMode',figure_mode,'FigureLayout',figure_layout, ...
        'OutputPath',output_path,'FigurePrefix', ...
        ['RepeatedSessions_BN' num2str(selected_networks(network_index))]);
end

%% 6) COLLECT THE STATISTICAL OUTPUTS

SESSION_STATS.SelectedNetworks = selected_networks;
SESSION_STATS.NumberConditions = number_conditions;
SESSION_STATS.NumberSessions = number_sessions;
SESSION_STATS.SessionPairs = session_pairs;
SESSION_STATS.ContrastLabels = session_contrast_labels;
SESSION_STATS.Time = time;
SESSION_STATS.PValues = p_values;
SESSION_STATS.AdjustedPValues = adjusted_p;
SESSION_STATS.Significant = significant;
SESSION_STATS.FDRCriticalP = critical_p;
SESSION_STATS.TStatistics = t_statistics;
SESSION_STATS.CohenDz = cohen_dz;
SESSION_STATS.SummaryTimeWindow = summary_time_window;
SESSION_STATS.WindowModels = window_models;
SESSION_STATS.WindowANOVA = window_anova;
SESSION_STATS.BehaviourModels = behaviour_models;
SESSION_STATS.Figures = statistical_figures;

disp('Repeated-session statistical example completed. Results are stored in SESSION_STATS.');
