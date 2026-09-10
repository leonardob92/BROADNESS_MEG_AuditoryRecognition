% ========================================================================
%  BROADNESS STATISTICS EXAMPLE 1: CONDITIONS WITHIN ONE PARTICIPANT GROUP
% ========================================================================
%
%  This script demonstrates participant-level comparisons between
%  experimental conditions using BROADNESS network time series.
%
%  The script loads one file per participant, concatenates the data along
%  the fourth dimension, and runs BROADNESS_NetworkEstimation. The resulting
%  network time series have dimensions:
%  time x networks x conditions x participants.
%
%  The participant is the independent unit of inference. Time-points are
%  tested in parallel and are never treated as independent participants.
%  Paired tests identify when specific pairs of conditions differ. An
%  optional analysis also illustrates how the network responses can be
%  related to one behavioural value measured for each participant.
%
%  Statistics and Machine Learning Toolbox is required for ttest and corr.
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
condition_pairs = nchoosek(1:number_conditions,2); %all pairwise condition contrasts
condition_labels = {'Memorized','NewT1','NewT2','NewT3','NewT4'};
alpha = 0.05;
minimum_significant_samples = 1;
if ~exist('behaviour', 'var')
    % Replace [] with one value per participant, ordered exactly as "list"
    % and therefore as the fourth dimension of "network_time_series".
    behaviour = [];
end
correlation_type = 'Pearson'; %'Pearson' or 'Spearman'

% To test this section without real behavioural data, uncomment the lines
% below. This creates a test variable from Brain Network 1, Condition 1.
% Use this only to check that the code works, not for scientific inference.
test_indices = time >= 0.30 & time <= 0.50;
behaviour = squeeze(mean(network_time_series(test_indices,1,1,:),1));
behaviour = behaviour(:);

figure_mode = 'show'; %'off', 'show', 'save', or 'both'
figure_layout = 'individual';

if length(condition_labels) ~= number_conditions
    error('Provide one entry in "condition_labels" for every condition.');
end
if ~ismember(lower(correlation_type), {'pearson','spearman'})
    error('"correlation_type" must be ''Pearson'' or ''Spearman''.');
end


%% 3) PAIRED CONDITION TESTS ACROSS TIME

run_condition_anova = 'on'; %'on' or 'off'

% One paired t-test is performed at every time-point, for every selected
% network and every pair of conditions. The tests are paired because the
% same participants contribute data to both conditions. These tests answer:
% "Between which conditions, and at which time-points, is there evidence of
% a difference?"
%
% This worked example uses paired t-tests. If the paired differences do not
% meet the intended parametric assumptions, signrank can be substituted at
% each time-point while retaining the same FDR and plotting workflow.

number_pairs = size(condition_pairs,1);
p_values = nan(number_timepoints,length(selected_networks),number_pairs);
t_statistics = nan(size(p_values));
cohen_dz = nan(size(p_values));

for network_index = 1:length(selected_networks)
    network = selected_networks(network_index);
    for pair = 1:number_pairs
        condition_1 = condition_pairs(pair,1);
        condition_2 = condition_pairs(pair,2);
        values_1 = reshape(network_time_series(:,network,condition_1,:), ...
            number_timepoints,number_participants);
        values_2 = reshape(network_time_series(:,network,condition_2,:), ...
            number_timepoints,number_participants);
        [~,p_values(:,network_index,pair),~,test_statistics] = ...
            ttest(values_1,values_2,'Dim',2);
        t_statistics(:,network_index,pair) = test_statistics.tstat;
        paired_differences = values_1 - values_2;
        cohen_dz(:,network_index,pair) = mean(paired_differences,2,'omitnan') ./ ...
            std(paired_differences,0,2,'omitnan');
        disp(pair)
    end
end
cohen_dz(~isfinite(cohen_dz)) = NaN;

% All selected time-points, networks, and condition contrasts form one FDR
% family here. Change the supplied family only with an a priori rationale.
[significant,adjusted_p,critical_p] = BROADNESS_FDRCorrection(p_values,alpha);

contrast_labels = cell(number_pairs,1);
for pair = 1:number_pairs
    contrast_labels{pair} = [condition_labels{condition_pairs(pair,1)} ...
        ' vs ' condition_labels{condition_pairs(pair,2)}];
end

% The optional repeated-measures ANOVA tests the overall condition effect
% at every time-point. With more than two conditions, Greenhouse-Geisser
% correction is applied before FDR correction across time and networks.
if strcmpi(run_condition_anova,'on')
    condition_anova = BROADNESS_RepeatedMeasuresANOVA( ...
        network_time_series(:,selected_networks,:,:),alpha);
elseif strcmpi(run_condition_anova,'off')
    condition_anova.Enabled = false;
    condition_anova.Reason = 'The optional condition ANOVA was switched off.';
    condition_anova.Alpha = alpha;
else
    error('"run_condition_anova" must be ''on'' or ''off''.');
end

% Collect the condition-comparison statistics from this section.
CONDITION_STATS.SelectedNetworks = selected_networks;
CONDITION_STATS.ConditionLabels = condition_labels;
CONDITION_STATS.Time = time;
CONDITION_STATS.Alpha = alpha;
CONDITION_STATS.PairedTests.ConditionPairs = condition_pairs;
CONDITION_STATS.PairedTests.ContrastLabels = contrast_labels;
CONDITION_STATS.PairedTests.PValues = p_values;
CONDITION_STATS.PairedTests.AdjustedPValues = adjusted_p;
CONDITION_STATS.PairedTests.Significant = significant;
CONDITION_STATS.PairedTests.FDRCriticalP = critical_p;
CONDITION_STATS.PairedTests.FDRMethod = 'Benjamini-Hochberg';
CONDITION_STATS.PairedTests.FDRFamily = ...
    'All selected time-points, networks, and condition pairs';
CONDITION_STATS.PairedTests.TStatistics = t_statistics;
CONDITION_STATS.PairedTests.CohenDz = cohen_dz;
CONDITION_STATS.ANOVA = condition_anova;

%% 4) PLOT THE FDR-SIGNIFICANT WINDOWS

statistical_figures = cell(length(selected_networks),1);
pairwise_colors = lines(max(1,number_pairs));
pairwise_codes = arrayfun(@num2str,1:number_pairs,'UniformOutput',false)';

if condition_anova.Enabled
    significance_colors = [0 0 0; pairwise_colors];
    significance_labels = [{'A'}; pairwise_codes];
    significance_key = table([{'A'}; pairwise_codes], ...
        [{'Overall condition effect (ANOVA)'}; contrast_labels], ...
        'VariableNames',{'Line','StatisticalEffect'});
else
    significance_colors = pairwise_colors;
    significance_labels = pairwise_codes;
    significance_key = table(pairwise_codes,contrast_labels, ...
        'VariableNames',{'Line','StatisticalEffect'});
end

disp('Significance-line key:')
disp(significance_key)

for network_index = 1:length(selected_networks)
    pairwise_windows = cell(number_pairs,1);
    for pair = 1:number_pairs
        pairwise_windows{pair} = BROADNESS_SignificantSamplesToWindows( ...
            significant(:,network_index,pair),time,minimum_significant_samples);
    end
    if condition_anova.Enabled
        anova_windows = BROADNESS_SignificantSamplesToWindows( ...
            condition_anova.Significant(:,network_index),time, ...
            minimum_significant_samples);
        significant_windows = [{anova_windows}; pairwise_windows];
    else
        significant_windows = pairwise_windows;
    end
    statistical_figures{network_index} = BROADNESS_Plot_ActivationTimeseries( ...
        network_time_series(:,selected_networks(network_index),:,:),time, ...
        'ConditionLabels',condition_labels, ...
        'SignificantWindows',significant_windows, ...
        'SignificanceColors',significance_colors, ...
        'SignificanceLabels',significance_labels, ...
        'SignificanceStyle','line', ...
        'SignificanceLineOffset',0.025, ...
        'FigureMode',figure_mode,'FigureLayout',figure_layout, ...
        'OutputPath',output_path,'FigurePrefix', ...
        ['WithinParticipants_BN' num2str(selected_networks(network_index))]);
end

CONDITION_STATS.SignificanceKey = significance_key;
CONDITION_STATS.Figures = statistical_figures;

%% 5) OPTIONAL ASSOCIATION WITH A BEHAVIOURAL MEASURE

% This optional analysis relates the network response to one behavioural
% value measured for each participant. At every time-point, the participants'
% network responses are correlated with their behavioural values. The
% calculation is performed separately for every selected network and every
% experimental condition. It therefore answers:
% "When, and in which network and condition, is network activity associated
% with individual differences in behaviour?"
%
% Pearson correlation tests a linear association. Spearman correlation can
% instead be selected for a rank-based monotonic association. All tested
% time-points, networks, and conditions are included in one FDR family.
% Leave "behaviour" empty to skip this analysis.

behaviour_correlations = nan(number_timepoints,length(selected_networks), ...
    number_conditions);
behaviour_p_values = nan(size(behaviour_correlations));
behaviour_adjusted_p = nan(size(behaviour_correlations));
behaviour_significant = false(size(behaviour_correlations));
behaviour_critical_p = NaN;
behaviour_correlation_figures = gobjects(0);
behaviour_correlation_files = {};

if ~isempty(behaviour)
    selected_time_series = network_time_series(:,selected_networks,:,:);
    [behaviour_correlations,behaviour_p_values,behaviour_adjusted_p, ...
        behaviour_significant,behaviour_critical_p] = ...
        BROADNESS_BehaviourCorrelation(selected_time_series,behaviour, ...
        alpha,correlation_type);

    % Plot the correlation coefficient across time. Thick portions of each
    % curve indicate samples that remain significant after FDR correction.
    if ~strcmpi(figure_mode,'off')
        behaviour_figure_settings = BROADNESS_FigureSettings(figure_mode, ...
            'individual',output_path,{'png'},'WithinParticipants', ...
            'BehaviourCorrelations');
        condition_colors = lines(number_conditions);
        for network_index = 1:length(selected_networks)
            fig = figure('Visible',behaviour_figure_settings.Visible, ...
                'Color','w','Position',[100 100 850 480]);
            ax = gca;
            hold(ax,'on')
            correlation_lines = gobjects(number_conditions,1);
            for condition = 1:number_conditions
                correlation_curve = behaviour_correlations(:,network_index,condition);
                correlation_lines(condition) = plot(ax,time,correlation_curve, ...
                    'Color',condition_colors(condition,:),'LineWidth',1.5);
                significant_curve = correlation_curve;
                significant_curve(~behaviour_significant(:,network_index,condition)) = NaN;
                plot(ax,time,significant_curve,'Color',condition_colors(condition,:), ...
                    'LineWidth',4,'HandleVisibility','off');
            end
            yline(ax,0,':','Color',[0.35 0.35 0.35], ...
                'HandleVisibility','off');
            xlabel(ax,'Time (s)')
            ylabel(ax,[correlation_type ' correlation'])
            title(ax,['Behaviour association - Brain Network ' ...
                num2str(selected_networks(network_index))])
            legend(ax,correlation_lines,condition_labels,'Location','best', ...
                'Box','off')
            xlim(ax,[time(1) time(end)])
            ylim(ax,[-1 1])
            grid(ax,'on')
            box(ax,'off')
            set(ax,'FontSize',12,'LineWidth',1,'TickDir','out')
            behaviour_correlation_figures(end+1,1) = fig; %#ok<SAGROW>
            saved_files = BROADNESS_FinalizeFigure(fig, ...
                behaviour_figure_settings,['BehaviourCorrelation_BN' ...
                num2str(selected_networks(network_index))]);
            behaviour_correlation_files = [behaviour_correlation_files; ...
                saved_files]; %#ok<AGROW>
        end

    end
end

% Collect the network-amplitude behavioural statistics from this section.
BEHAVIOUR_STATS.Enabled = ~isempty(behaviour);
BEHAVIOUR_STATS.SelectedNetworks = selected_networks;
BEHAVIOUR_STATS.ConditionLabels = condition_labels;
BEHAVIOUR_STATS.Time = time;
BEHAVIOUR_STATS.Behaviour = behaviour;
BEHAVIOUR_STATS.CorrelationType = correlation_type;
BEHAVIOUR_STATS.Alpha = alpha;
BEHAVIOUR_STATS.Correlations = behaviour_correlations;
BEHAVIOUR_STATS.PValues = behaviour_p_values;
BEHAVIOUR_STATS.AdjustedPValues = behaviour_adjusted_p;
BEHAVIOUR_STATS.Significant = behaviour_significant;
BEHAVIOUR_STATS.FDRCriticalP = behaviour_critical_p;
BEHAVIOUR_STATS.FDRMethod = 'Benjamini-Hochberg';
BEHAVIOUR_STATS.FDRFamily = ...
    'All selected time-points, networks, and conditions';
BEHAVIOUR_STATS.Figures = behaviour_correlation_figures;
BEHAVIOUR_STATS.FigureFiles = behaviour_correlation_files;

%% 6) OPTIONAL BEHAVIOURAL ASSOCIATION WITH PHASE-SPACE OR RQA METRICS

%%% ------------------- USER SETTINGS ------------------- %%%

phase_behaviour_analysis = 'RQA'; %'PhaseSpace' or 'RQA'
selected_phase_metrics = {'Dispersion','Speed'}; %used only for 'PhaseSpace'

% With 'PhaseSpace', select one or more time-resolved metrics from:
% 'Dispersion', 'Speed', 'Acceleration', 'BaselineDisplacement', or
% 'CumulativePathLength'. With 'RQA', all eight RQA metrics are analysed.
%
% The correlations are calculated separately for every condition. FDR is
% applied jointly across the selected metrics, conditions, and time-points
% (phase-space analysis), or across all RQA metrics and conditions.

PHASE_RQA_BEHAVIOUR.Enabled = false;
PHASE_RQA_BEHAVIOUR.Behaviour = behaviour;
PHASE_RQA_BEHAVIOUR.CorrelationType = correlation_type;
PHASE_RQA_BEHAVIOUR.Alpha = alpha;
PHASE_RQA_BEHAVIOUR.FDRMethod = 'Benjamini-Hochberg';

if ~isempty(behaviour)
    if ~ismember(lower(phase_behaviour_analysis), {'phasespace','rqa'})
        error('"phase_behaviour_analysis" must be ''PhaseSpace'' or ''RQA''.');
    end

    % Compute the participant-level phase space and RQA only once. If this
    % output already exists in the workspace, it is reused.
    if ~exist('RQA_BROADNESS','var')
        RQA_BROADNESS = BROADNESS_PhaseSpace_RQA(BROADNESS, ...
            'principalcomps',selected_networks,'figuremode','off');
    end

    PHASE_RQA_BEHAVIOUR.Enabled = true;
    PHASE_RQA_BEHAVIOUR.Analysis = phase_behaviour_analysis;
    PHASE_RQA_BEHAVIOUR.Conditions = condition_labels;

    if strcmpi(phase_behaviour_analysis,'PhaseSpace')
        available_phase_metrics = {'Dispersion','Speed','Acceleration', ...
            'BaselineDisplacement','CumulativePathLength'};
        if ~iscell(selected_phase_metrics) || isempty(selected_phase_metrics) || ...
                any(~ismember(selected_phase_metrics,available_phase_metrics))
            error(['"selected_phase_metrics" must contain one or more valid ' ...
                'time-resolved phase-space metrics.']);
        end

        % Obtain the participant-level phase-space metrics without running
        % the inferential tests in BROADNESS_PhaseSpaceStatistics.
        PHASE_STATS = BROADNESS_PhaseSpaceStatistics(RQA_BROADNESS, ...
            'conditionnames',condition_labels,'statistics','off');
        phase_time = PHASE_STATS.Time(:);
        number_phase_metrics = length(selected_phase_metrics);
        phase_metric_data = nan(length(phase_time),number_phase_metrics, ...
            number_conditions,number_participants);

        for metric_index = 1:number_phase_metrics
            phase_metric_data(:,metric_index,:,:) = ...
                PHASE_STATS.Metrics.(selected_phase_metrics{metric_index});
        end

        [correlations,p_values_phase,adjusted_p_phase,significant_phase, ...
            critical_p_phase] = BROADNESS_BehaviourCorrelation( ...
            phase_metric_data,behaviour,alpha,correlation_type);

        PHASE_RQA_BEHAVIOUR.MetricNames = selected_phase_metrics;
        PHASE_RQA_BEHAVIOUR.Time = phase_time;
        PHASE_RQA_BEHAVIOUR.FDRFamily = ...
            'Selected phase-space metrics, conditions, and time-points';

    else
        % Each participant has one RQA table containing all conditions and
        % all eight metrics. These summary metrics do not have a time axis.
        rqa_metric_names = RQA_BROADNESS.RQA_metrics{1}.Properties.VariableNames;
        number_rqa_metrics = length(rqa_metric_names);
        rqa_metric_data = nan(number_rqa_metrics,number_conditions, ...
            number_participants);

        for participant = 1:number_participants
            rqa_metric_data(:,:,participant) = ...
                table2array(RQA_BROADNESS.RQA_metrics{participant})';
        end

        [correlations,p_values_phase,adjusted_p_phase,significant_phase, ...
            critical_p_phase] = BROADNESS_BehaviourCorrelation( ...
            rqa_metric_data,behaviour,alpha,correlation_type);

        PHASE_RQA_BEHAVIOUR.MetricNames = rqa_metric_names;
        PHASE_RQA_BEHAVIOUR.Time = [];
        PHASE_RQA_BEHAVIOUR.FDRFamily = 'All RQA metrics and conditions';
    end

    PHASE_RQA_BEHAVIOUR.Correlations = correlations;
    PHASE_RQA_BEHAVIOUR.PValues = p_values_phase;
    PHASE_RQA_BEHAVIOUR.AdjustedPValues = adjusted_p_phase;
    PHASE_RQA_BEHAVIOUR.Significant = significant_phase;
    PHASE_RQA_BEHAVIOUR.FDRCriticalP = critical_p_phase;
end

disp(['Within-participant statistical example completed. Results are stored ' ...
    'in CONDITION_STATS, BEHAVIOUR_STATS, and PHASE_RQA_BEHAVIOUR.']);

%%
