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
%  HOW TO USE THIS EXAMPLE:
%  1. Set the paths and participant-file pattern in Sections 0-1.
%  2. Choose networks, conditions and statistics in Section 2.
%  3. Set figure options in Section 3; spatial options in Section 4.
%  4. Optionally supply independent behavioural scores in Sections 5-6.
%  5. Inspect CONDITION_STATS, SPATIAL_STATS, BEHAVIOUR_STATS and
%     PHASE_RQA_BEHAVIOUR. Participant order must match the loaded files.
%
%  Requires the BROADNESS toolbox and MATLAB Statistics Toolbox.
% ========================================================================


%% 0) STARTUP

% Set the paths below, then run the script from top to bottom.
% Each analysis section begins with the settings relevant to that section.

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

% Load one file per participant; each file must contain "Data".
% Input dimensions: sources x time x conditions x participants.
list = dir(fullfile(data_path, 'SUBJ*.mat'));
data = [];
for subi = 1:length(list) %over participants
    load(fullfile(list(subi).folder, list(subi).name)) %loading data for each participant
    data = cat(4,data,Data); %concatenating data
    disp(subi)
end
data(:,777:end,:,:) = []; % example-specific truncation; adjust for your dataset
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


%% 2) PAIRED CONDITION TESTS ACROSS TIME

% USER SETTINGS: select the networks, conditions and correction method.
% Condition pairs are generated automatically (e.g. 5 conditions = 10 pairs).
selected_networks = [1 2];
condition_labels = {'Memorized','NewT1','NewT2','NewT3','NewT4'};
condition_pairs = nchoosek(1:number_conditions,2);
alpha = 0.05;
correction_method = 'cluster'; % 'FDR' or 'cluster'
cluster_permutations = 5000; % used only for cluster correction
cluster_alpha = 0.05; % cluster-forming threshold
run_condition_anova = 'on'; % 'on' or 'off'; overall condition test

if length(condition_labels) ~= number_conditions
    error('Provide one entry in "condition_labels" for every condition.');
end

% HOW TO READ RESULTS: CONDITION_STATS.PairedTests contains the paired
% comparisons (including CohenDz); CONDITION_STATS.ANOVA contains the
% overall condition test. Cluster p-values refer to whole clusters.


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

% Correction for multiple comparisons
if strcmpi(correction_method,'FDR')
    [significant,adjusted_p,critical_p] = ...
        BROADNESS_FDRCorrection(p_values,alpha);
    cluster_stats = [];
    correction_name = 'Benjamini-Hochberg';
    correction_family = ...
        'All selected time-points, networks, and condition pairs';
elseif strcmpi(correction_method,'cluster')
    [significant,adjusted_p,critical_p,cluster_stats] = ...
        BROADNESS_TemporalClusterCorrection( ...
        network_time_series(:,selected_networks,:,:), ...
        condition_pairs,alpha, ...
        'permutations',cluster_permutations, ...
        'clusteralpha',cluster_alpha, ...
        'time',time);
    correction_name = 'Temporal cluster-mass permutation';
    correction_family = ...
        'All time-points, separately per network and condition pair';
else
    error('correction_method must be ''FDR'' or ''cluster''.');
end
contrast_labels = cell(number_pairs,1);
for pair = 1:number_pairs
    contrast_labels{pair} = [condition_labels{condition_pairs(pair,1)} ...
        ' vs ' condition_labels{condition_pairs(pair,2)}];
end

% The optional repeated-measures ANOVA tests the overall condition effect
% at every time-point. The FDR option uses Greenhouse-Geisser-adjusted
% p-values when appropriate; the cluster option uses condition permutations.
if strcmpi(run_condition_anova,'on')
    if strcmpi(correction_method,'FDR')
        condition_anova = BROADNESS_RepeatedMeasuresANOVA( ...
            network_time_series(:,selected_networks,:,:),alpha);
    elseif strcmpi(correction_method,'cluster')
        [anova_significant,anova_corrected_p,~,anova_cluster_stats] = ...
            BROADNESS_TemporalClusterCorrection_ANOVA( ...
            network_time_series(:,selected_networks,:,:),alpha, ...
            'permutations',cluster_permutations, ...
            'clusteralpha',cluster_alpha, ...
            'time',time);
        condition_anova.Enabled = true;
        condition_anova.Significant = anova_significant;
        condition_anova.FStatistics = anova_cluster_stats.FStatistics;
        condition_anova.ClusterCorrectedP = anova_corrected_p;
        condition_anova.ClusterResults = anova_cluster_stats;
        condition_anova.CorrectionMethod = 'Temporal cluster-mass permutation';
    else
        error('Unknown correction_method.');
    end
elseif strcmpi(run_condition_anova,'off')
    condition_anova.Enabled = false;
    condition_anova.Reason = 'The optional condition ANOVA was switched off.';
    condition_anova.Alpha = alpha;
else
    error('"run_condition_anova" must be ''on'' or ''off''.');
end

% Collect the condition-comparison statistics from this section.
% General analysis settings
CONDITION_STATS.Settings.SelectedNetworks = selected_networks;
CONDITION_STATS.Settings.ConditionLabels = condition_labels;
CONDITION_STATS.Settings.Time = time;
CONDITION_STATS.Settings.Alpha = alpha;
CONDITION_STATS.Settings.RunConditionANOVA = run_condition_anova;
CONDITION_STATS.Settings.CorrectionMethod = correction_method;
CONDITION_STATS.Settings.ClusterPermutations = cluster_permutations;
CONDITION_STATS.Settings.ClusterAlpha = cluster_alpha;
% Paired condition comparisons
CONDITION_STATS.PairedTests.ConditionPairs = condition_pairs;
CONDITION_STATS.PairedTests.ContrastLabels = contrast_labels;
CONDITION_STATS.PairedTests.PValues = p_values;
CONDITION_STATS.PairedTests.AdjustedPValues = adjusted_p;
CONDITION_STATS.PairedTests.Significant = significant;
CONDITION_STATS.PairedTests.FDRCriticalP = critical_p;
CONDITION_STATS.PairedTests.CorrectionMethod = correction_name;
CONDITION_STATS.PairedTests.CorrectionFamily = correction_family;
CONDITION_STATS.PairedTests.ClusterResults = cluster_stats;
CONDITION_STATS.PairedTests.TStatistics = t_statistics;
CONDITION_STATS.PairedTests.CohenDz = cohen_dz;
% Overall condition ANOVA
CONDITION_STATS.ANOVA = condition_anova;

%% 3) PLOT THE SIGNIFICANT TIME WINDOWS

% USER SETTINGS: choose how figures are displayed and how short intervals
% may be. One sample is enough to display an interval by default.
figure_mode = 'show'; % 'off', 'show', 'save', or 'both'
figure_layout = 'individual';
minimum_significant_samples = 1;

% Curves show mean network responses by condition. The 'A' line marks the
% overall ANOVA; numbered lines correspond to condition pairs in the key.
% FDR marks individual samples; cluster correction tests whole intervals.

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

% Save the figure settings alongside the corresponding figure results.
CONDITION_STATS.Settings.MinimumSignificantSamples = minimum_significant_samples;
CONDITION_STATS.Settings.FigureMode = figure_mode;
CONDITION_STATS.Settings.FigureLayout = figure_layout;

% Collect figures and their significance-line descriptions
CONDITION_STATS.Figures.SignificanceKey = significance_key;
CONDITION_STATS.Figures.ConditionComparisons = statistical_figures;

%% 4) CONDITION DIFFERENCES IN SPATIAL ACTIVATION PATTERNS

% Participant spatial activation patterns can be compared between
% conditions at every source. The function uses unthresholded maps and
% cluster-mass permutation correction. It also calculates map strength,
% weighted centroid, spatial dispersion, and hemispheric lateralization.
%
% This example tests Network 1, Conditions 1 versus 2, using 100
% permutations for a quick check; use at least 5000 for final inference.
% To test all condition pairs, set "spatial_condition_pairs" to
% nchoosek(1:number_conditions,2).

% USER SETTINGS: these choices apply only to the spatial analysis.
% This section uses participant spatial maps, not the temporal t-tests.
selected_spatial_networks = 1;
spatial_condition_pairs = [1 2];
spatial_permutations = 100;
spatial_map_normalization = 'none'; %'none' or 'rms'
run_spatial_omnibus = 'on'; %'on' tests the overall condition effect

load(fullfile(project_path,'BROADNESS_External', ...
    'MNI152_8mm_coord_dyi.mat'),'MNI8');

[SPATIAL_STATS,SPATIAL_FIGURES] = BROADNESS_SpatialPatternStatistics( ...
    BROADNESS,'principalcomps',selected_spatial_networks, ...
    'conditionpairs',spatial_condition_pairs, ...
    'conditionnames',condition_labels,'mni_coords',MNI8, ...
    'mapnormalization',spatial_map_normalization, ...
    'permutations',spatial_permutations,'omnibus',run_spatial_omnibus, ...
    'figuremode',figure_mode,'figurelayout','both', ...
    'OutputPath',output_path,'figureprefix','WithinParticipants', ...
    'clusteralpha', 0.05, 'alpha', 0.05);


%% 5) OPTIONAL ASSOCIATION WITH A BEHAVIOURAL MEASURE

% USER SETTINGS: one external behavioural score per participant, in the
% same order as "list". Leave [] to skip this section and Section 6.
behaviour = [];
correlation_type = 'Pearson'; % 'Pearson' or 'Spearman'

% DEMONSTRATION ONLY: Uncomment to check the workflow with artificial
% behaviour derived from the MEG data. NEVER use it for scientific inference.
% test_indices = time >= 0.30 & time <= 0.50;
% behaviour = squeeze(mean(network_time_series(test_indices,1,1,:),1));
% behaviour = behaviour(:);

if ~ismember(lower(correlation_type), {'pearson','spearman'})
    error('"correlation_type" must be ''Pearson'' or ''Spearman''.');
end
if ~isempty(behaviour) && numel(behaviour) ~= number_participants
    error('Provide exactly one behavioural value per participant.');
end

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

% USER SETTINGS: choose phase-space metrics or RQA summaries.
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
    'in CONDITION_STATS, SPATIAL_STATS, BEHAVIOUR_STATS, and ' ...
    'PHASE_RQA_BEHAVIOUR.']);

%%
