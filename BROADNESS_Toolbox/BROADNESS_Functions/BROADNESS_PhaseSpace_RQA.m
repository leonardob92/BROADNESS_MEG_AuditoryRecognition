function [RQA_BROADNESS, FIGURES] = BROADNESS_PhaseSpace_RQA(BROADNESS, varargin)
%%
% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  RECURRENCE QUANTIFICATION ANALYSIS (RQA)
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
%  This function computes RQA metrics on BROADNESS-derived time series of
%  brain networks. It reconstructs the system’s phase space using selected
%  principal components, computes recurrence plots (with and without thresholding),
%  and extracts standard Recurrence Quantification Analysis (RQA) measures describing the system's temporal dynamics.
%  Here, phase space is the multivariate representation in which each point
%  describes the simultaneous activity of the selected brain networks at
%  one time point.
%
%  Specifically, it:
%   - Projects BROADNESS time series (brain networks) into a phase space using specified 
%     principal components (PCs). Each dimension in the phase space
%     correspond to a principal component.
%   - Optionally displays 2D/3D animated scatter plots of phase space
%   - Computes recurrence plots (multivariate distance matrices) and thresholded 
%     recurrence plots based on a user-defined threshold
%   - Extracts standard RQA metrics:
%       * Recurrence Rate (RR)
%       * Mean diagonal length (L)
%       * Determinism (DET)
%       * Entropy (ENTR)
%       * Trapping Time (TT)
%       * Laminarity (LAM)
%       * Maximal vertical line length (V_max)
%       * Divergence (DIV)
%   - Returns results in a structured output containing plots, matrices, 
%     and the metrics table
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - BROADNESS                          : Structure from BROADNESS_NetworkEstimation
%      - .TimeSeries_BrainNetworks      : 2D, 3D or 4D matrix (time × components × [conditions] x [participants])
%      - .Time                          : Vector of time points (seconds)
%
%  - Optional arguments (name-value pairs):
%      - 'principalcomps'               : Vector of PC indices to use in phase space (default: 1:10)
%      - 'timeinterval'                 : [start_time end_time] in seconds for analysis (default: full range)
%      - 'threshold'                    : Fraction of max distance to define recurrences (default: 0.1)
%      - 'normalization'                : Scaling applied to the selected phase-space dimensions:
%                                         'none' (default) preserves the original PCA/ICA scores;
%                                         'pooled_zscore' uses one mean and standard deviation per
%                                         selected component, pooled across the analysed time-points,
%                                         conditions, and participants
%      - 'theiler_window'               : Theiler window in samples. Use 0 to exclude only the main diagonal,
%                                         or a positive integer to exclude a wider diagonal band.
%                                         Default: [] (disabled; main diagonal retained)
%      - 'video'                        : 'on' or 'off' to show animated phase space plot (default: 'off')
%      - 'figure'                       : 'on' or 'off' to show the figures (default: 'off')
%      - 'figuremode'                   : 'off', 'show', 'save', or 'both'. If omitted,
%                                         the legacy 'figure' option is used (default: 'off')
%      - 'figurelayout'                 : 'individual', 'summary', or 'both' (default: 'individual')
%      - 'OutputPath'                   : Base output folder required when figures are saved
%      - 'figureformats'                : 'png', 'pdf', 'fig', or a cell array (default: {'png'})
%      - 'figureprefix'                 : Optional prefix for saved figure filenames
%      - 'outpath'                      : Deprecated alias for 'OutputPath', retained for
%                                         compatibility with previous BROADNESS scripts
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - RQA_BROADNESS                              : Structure with RQA results
%      - .PhaseSpace.Time                      : Analysed time points in seconds
%      - .PhaseSpace.PCs                       : Brain-network components used as dimensions
%      - .PhaseSpace.ParticipantCoordinates    : Cell array of participant trajectories
%      - .PhaseSpace.MeanCoordinates           : Cell array of participant-averaged trajectories
%      - .PhaseSpace.Normalization              : Applied method and common center/scale parameters
%      - .RecurrencePlots.DistMat               : Cell array of distance matrices
%      - .RecurrencePlots.RecurPlot             : Cell array of recurrence plots (i.e., thresholded distance matrices)
%      - .RQA_metrics                           : Table of 8 RQA measures
%  - FIGURES                                    : Visible figure handles and saved figure paths
%
%      Please, note that this output will be generated for each experimental condition and participant,
%      if the data was originally provided in such format. 
%
% ------------------------------------------------------------------------
%
% ========================================================================
%  AUTHORS:
%  Chiara Malvaso, Mattia Rosso & Leonardo Bonetti 
%  chiara.malvaso@studio.unibo.it
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Department of Physics, University of Bologna
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 11/10/2025
%
% ========================================================================
%
%
% NOTE: This function uses code from the Cross Recurrence Plot Toolbox by Norbert Marwan:
% Norbert Marwan (2025). Cross Recurrence Plot Toolbox for MATLAB
% (https://www.mathworks.com/matlabcentral/fileexchange/6170-cross-recurrence-plot-toolbox)
% MATLAB Central File Exchange. Retrieved July 05, 2025.
%
%%








%% ----------------------------- Parse inputs -----------------------------

disp('Checking inputs')

% Defaults
opts = struct('principalcomps', 1:2, 'timeinterval', [], 'threshold', 0.1, ...
    'normalization', 'none', 'theiler_window', [], 'video', 'off', 'figure', 'off', ...
    'figuremode', [], 'figurelayout', 'individual', 'outputpath', [], ...
    'outpath', [], ...
    'figureformats', {{'png'}}, 'figureprefix', '');
opts = parse_name_value_pairs(opts, varargin{:});

% Assign to readable internal names
PCs          = opts.principalcomps;
time_seconds = opts.timeinterval;
eps          = opts.threshold;
normalization = lower(char(string(opts.normalization)));
theiler_window = opts.theiler_window;
video        = opts.video;
figurel      = opts.figure;

if isempty(opts.figuremode)
    if strcmpi(figurel, 'on')
        opts.figuremode = 'show';
    else
        opts.figuremode = 'off';
    end
end
outputPath = opts.outputpath;
if isempty(outputPath)
    outputPath = opts.outpath;
elseif ~isempty(opts.outpath) && ...
        ~strcmp(char(string(opts.outpath)), char(string(outputPath)))
    warning(['The deprecated ''outpath'' value is ignored when ' ...
        '''OutputPath'' is provided.']);
end
figureSettings = BROADNESS_FigureSettings(opts.figuremode, opts.figurelayout, ...
    outputPath, opts.figureformats, opts.figureprefix, 'PhaseSpace_RQA');
figureHandles = gobjects(0);
figureFiles = {};

% --- Required fields FIRST (so 'time' is defined before we use it) ---
if isfield(BROADNESS, 'TimeSeries_BrainNetworks')
    TimeSeries = BROADNESS.TimeSeries_BrainNetworks;
else
    error('Invalid input structure: field "TimeSeries_BrainNetworks" is required.');
end
if isfield(BROADNESS, 'Time')
    time = BROADNESS.Time;
else
    error('Invalid input structure: field "Time" is required.');
end
if size(TimeSeries,1) ~= numel(time)
    error('The first dimension of "TimeSeries_BrainNetworks" must match the length of "Time".');
end

% --- Validate optional args ---
if ~(isnumeric(PCs) && isvector(PCs))
    error('"principalcomps" must be a numeric vector.');
end
if ~isnumeric(eps) || ~isscalar(eps) || ~isfinite(eps)
    error('"threshold" must be a finite numeric scalar.');
end
if ~ismember(normalization, {'none','pooled_zscore'})
    error('"normalization" must be ''none'' or ''pooled_zscore''.');
end
if ~isempty(theiler_window) && (~isnumeric(theiler_window) || ~isscalar(theiler_window) || ...
        ~isfinite(theiler_window) || theiler_window < 0 || fix(theiler_window) ~= theiler_window)
    error('"theiler_window" must be empty or a nonnegative integer number of samples.');
end
if ~(ischar(video) || isstring(video)) || ~ismember(lower(string(video)), ["on","off"])
    error('"video" must be ''on'' or ''off''.');
end
if ~(ischar(figurel) || isstring(figurel)) || ~ismember(lower(string(figurel)), ["on","off"])
    error('"figure" must be ''on'' or ''off''.');
end

% ---- Validate/normalize time interval ----
if isempty(time_seconds)
    time_seconds = [time(1) time(end)];  % default full span
else
    if ~isnumeric(time_seconds) || numel(time_seconds) ~= 2 || any(~isfinite(time_seconds))
        error('"timeinterval" must be a numeric 1x2 vector: [start end] in seconds.');
    end
    if time_seconds(1) > time_seconds(2)
        time_seconds = time_seconds([2 1]); % be forgiving
    end
    % Clip to available range
    time_seconds(1) = max(time_seconds(1), time(1));
    time_seconds(2) = min(time_seconds(2), time(end));
    if time_seconds(1) >= time_seconds(2)
        error('Requested "timeinterval" is outside available data range.');
    end
end

% Indices for the chosen window (start >=, end <=). The small tolerance
% retains endpoints that differ only because of floating-point precision.
time_tolerance = max(10*builtin('eps',max(1,max(abs(time)))), median(diff(time))*1e-9);
timemin = find(time >= time_seconds(1)-time_tolerance, 1, 'first');
timemax = find(time <= time_seconds(2)+time_tolerance, 1, 'last');
reduced_time_idx = timemin:timemax;

if ~isempty(theiler_window)
    if theiler_window >= length(reduced_time_idx)-1
        error('"theiler_window" must leave at least one pair of eligible time-points.');
    end
    time_indices = 1:length(reduced_time_idx);
    theiler_mask = abs(time_indices.'-time_indices) > theiler_window;
end

if any(PCs < 1) || any(PCs > size(TimeSeries,2)) || any(fix(PCs) ~= PCs)
    error('"principalcomps" contains indices outside the available brain networks.');
end

%% ------------------- Compute phase space coordinates --------------------
% Now supports 4th dim = participants.
% Video (and downstream code using `phase_space`) uses the participant-AVERAGED TimeSeries.

% If TimeSeries is 3D, promote to 4D with singleton participants
if ndims(TimeSeries) == 3
    TimeSeries = reshape(TimeSeries, size(TimeSeries,1), size(TimeSeries,2), size(TimeSeries,3), 1);
end

nT   = size(TimeSeries,1);
nPC  = size(TimeSeries,2);
nCond = size(TimeSeries,3);
nPart = size(TimeSeries,4);

% Optionally place all selected dimensions on a common standardized scale.
% Crucially, the same center and scale are used for every condition and
% participant, so between-condition and between-group differences are not
% removed by condition-specific or participant-specific normalization.
normalizationCenter = zeros(1,length(PCs));
normalizationScale = ones(1,length(PCs));
if strcmp(normalization, 'pooled_zscore')
    for pcColumn = 1:length(PCs)
        pooledValues = TimeSeries(reduced_time_idx, PCs(pcColumn), :, :);
        pooledValues = pooledValues(:);
        if any(~isfinite(pooledValues))
            error(['Pooled z-score normalization requires finite values in ' ...
                'all selected time series.']);
        end
        normalizationCenter(pcColumn) = mean(pooledValues);
        normalizationScale(pcColumn) = std(pooledValues, 0);
        if ~isfinite(normalizationScale(pcColumn)) || ...
                normalizationScale(pcColumn) <= builtin('eps', ...
                max(1,max(abs(pooledValues))))
            error(['Pooled z-score normalization cannot be applied because ' ...
                'brain network ' num2str(PCs(pcColumn)) ...
                ' has zero or near-zero variance.']);
        end
        TimeSeries(:,PCs(pcColumn),:,:) = ...
            (TimeSeries(:,PCs(pcColumn),:,:) - normalizationCenter(pcColumn)) ...
            ./ normalizationScale(pcColumn);
    end
end

% 1) Per-participant phase spaces (so nothing is lost if you need them later)
phase_space_participants = cell(nCond, nPart);
for part = 1:nPart            % over participants
    for cond = 1:nCond        % over conditions
        % initialize the matrix that will contain the coordinates in the phase space
        phase_space_temp = zeros(length(reduced_time_idx), length(PCs));
        for t = 1:length(reduced_time_idx)   % over time points
            for cc = 1:length(PCs)           % over components
                phase_space_temp(t,cc) = TimeSeries(reduced_time_idx(t), PCs(cc), cond, part);
            end
        end
        phase_space_participants{cond, part} = phase_space_temp;
    end
end

% 2) Participant-AVERAGED TimeSeries (required for video display per your note)
TimeSeries_avg = mean(TimeSeries, 4); % average across participants -> 3D: [time x PC x cond]

% 3) Averaged phase_space (keep original variable name & shape so downstream code is unchanged)
phase_space = cell(nCond,1);
for cond = 1:nCond % over conditions
    phase_space_temp = zeros(length(reduced_time_idx), length(PCs));
    for t = 1:length(reduced_time_idx)      % over time points
        for cc = 1:length(PCs)              % over components
            phase_space_temp(t,cc) = TimeSeries_avg(reduced_time_idx(t), PCs(cc), cond);
        end
    end
    phase_space{cond} = phase_space_temp;
end

%% ---------------------- Display phase space video -----------------------

% will now use the averaged `phase_space`
if strcmp(video, 'on')
    if length(PCs) == 2
        targetTicks = 10; % how many ticks you want on the colorbar
        timeRange = timemax - timemin;
        X = round(timeRange / (targetTicks - 1));
        scatsize = 30;
        for cond = 1:size(phase_space,1) % over conditions
            figure
            xlim([min(phase_space{cond}(:,1)) max(phase_space{cond}(:,1))])
            ylim([min(phase_space{cond}(:,2)) max(phase_space{cond}(:,2))])
            set(gcf,'color','w')
            % legend('show')
            grid minor
            box on
            c = colorbar;
            colormap('jet')
            bumba = time(timemin:X:timemax);
            c.Ticks = linspace(0,1,length(bumba));
            c.TickLabels = bumba;
            xlabel('Brain network 1')
            ylabel('Brain network 2')
            nscaz = length(phase_space{cond}(:,1));
            cmap = jet(nscaz);
            title(['Phase space plot - Condition ' num2str(cond)])
            for ii = 1:nscaz
                hold on
                scatter(phase_space{cond}(ii,1),phase_space{cond}(ii,2),scatsize,cmap(ii,:),'filled')
                pause(0.01)
            end
        end
    else
        if length(PCs) == 3
            targetTicks = 10; % how many ticks you want on the colorbar
            timeRange = timemax - timemin;
            X = round(timeRange / (targetTicks - 1));
            scatsize = 30;
            for cond = 1:size(phase_space,1) % over conditions
                figure
                view(3)          % Default 3D view
                axis vis3d       % Keep aspect ratio fixed during rotation

                xlim([min(phase_space{cond}(:,1)) max(phase_space{cond}(:,1))])
                ylim([min(phase_space{cond}(:,2)) max(phase_space{cond}(:,2))])
                zlim([min(phase_space{cond}(:,3)) max(phase_space{cond}(:,3))])
                set(gcf,'color','w')
                % legend('show')
                grid minor
                box on
                c = colorbar;
                colormap('jet')
                bumba = time(timemin:X:timemax);
                c.Ticks = linspace(0,1,length(bumba));
                c.TickLabels = bumba;
                xlabel('Brain network 1')
                ylabel('Brain network 2')
                zlabel('Brain network 3')
                nscaz = length(phase_space{cond}(:,1));
                cmap = jet(nscaz);
                title(['Phase space plot - Condition ' num2str(cond)])
                for ii = 1:nscaz
                    hold on
                    scatter3(phase_space{cond}(ii,1),phase_space{cond}(ii,2),phase_space{cond}(ii,3), scatsize, 'MarkerFaceColor', cmap(ii,:), 'MarkerEdgeColor', cmap(ii,:))
                    pause(0.01)
                end
            end

        else
            warning('Phase space video is available only in 2 or 3 dimensions. ')
        end
    end
end

%% ----------------------- Compute Recurrence Plot ------------------------
% Computed from averaged `phase_space` so it matches the video

% Initializing the cell array that will contain the distance matrix for each condition 
DM = cell(size(phase_space,1), 1);

% Initializing the cell array that will contain the thresholded recurrence plot for each condition 
RP = cell(size(phase_space,1), 1);

for cc = 1:size(phase_space,1) %over conditions
   
   disp(['Computing recurrence plot for condition ' num2str(cc)])
   RP_temp = zeros(size(phase_space{cc},1));
   RP_thresh_temp = zeros(size(RP_temp,1));
    for ii = 1:size(phase_space{cc},1) %over time-points
        for jj = 1:size(phase_space{cc},1) %over time-points
            % compute the distance between every couple of points in the phase space
            RP_temp(ii,jj) = norm(phase_space{cc}(ii,:) - phase_space{cc}(jj,:)); %distance between the n components for time-points ii and jj
        end
    end
   RP_thresh_temp(RP_temp<max(RP_temp(:))*eps) = 1; %recurrent values
   if ~isempty(theiler_window)
       RP_thresh_temp(~theiler_mask) = 0; %exclude temporally adjacent recurrence points
   end
   DM{cc} = RP_temp;
   RP{cc} = RP_thresh_temp;
%    clear RP_temp
%    clear RP_thresh_temp
end
    
%% ------------------------- Recurrence plots -----------------------------

if ~strcmp(figureSettings.Mode, 'off')
    for cond = 1:numel(phase_space)
        conditionName = ['Condition_' num2str(cond,'%02d')];

        if figureSettings.MakeIndividual
            fig = figure('Visible', figureSettings.Visible, 'Color', 'w');
            plot_phase_space(gca, phase_space{cond}, time(reduced_time_idx), cond);
            figureFiles = [figureFiles; BROADNESS_FinalizeFigure(fig, figureSettings, ...
                [conditionName '_PhaseSpace'])]; %#ok<AGROW>
            if figureSettings.Show, figureHandles(end+1) = fig; end %#ok<AGROW>

            fig = figure('Visible', figureSettings.Visible, 'Color', 'w');
            plot_recurrence_matrix(gca, DM{cond}, time(reduced_time_idx), ...
                ['Condition ' num2str(cond) ' — Distance matrix'], false);
            figureFiles = [figureFiles; BROADNESS_FinalizeFigure(fig, figureSettings, ...
                [conditionName '_DistanceMatrix'])]; %#ok<AGROW>
            if figureSettings.Show, figureHandles(end+1) = fig; end %#ok<AGROW>

            fig = figure('Visible', figureSettings.Visible, 'Color', 'w');
            plot_recurrence_matrix(gca, RP{cond}, time(reduced_time_idx), ...
                ['Condition ' num2str(cond) ' — Recurrence plot (' num2str(eps*100) '%)'], true);
            figureFiles = [figureFiles; BROADNESS_FinalizeFigure(fig, figureSettings, ...
                [conditionName '_RecurrenceThresholded'])]; %#ok<AGROW>
            if figureSettings.Show, figureHandles(end+1) = fig; end %#ok<AGROW>
        end

        if figureSettings.MakeSummary
            fig = figure('Visible', figureSettings.Visible, 'Color', 'w', ...
                'Position', [100 100 1350 430]);
            layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
            title(layout, ['Phase-Space RQA — Condition ' num2str(cond)]);
            plot_phase_space(nexttile(layout), phase_space{cond}, time(reduced_time_idx), cond);
            plot_recurrence_matrix(nexttile(layout), DM{cond}, time(reduced_time_idx), ...
                'Distance matrix', false);
            plot_recurrence_matrix(nexttile(layout), RP{cond}, time(reduced_time_idx), ...
                ['Recurrence plot (' num2str(eps*100) '%)'], true);
            figureFiles = [figureFiles; BROADNESS_FinalizeFigure(fig, figureSettings, ...
                [conditionName '_Summary'])]; %#ok<AGROW>
            if figureSettings.Show, figureHandles(end+1) = fig; end %#ok<AGROW>
        end
    end
end

%% --------------------- PER-PARTICIPANT RP + METRICS ---------------------

% New block: compute recurrence plots and RQA metrics for each participant,
% using the per-participant phase spaces we already computed earlier.

% Containers
RP_participants        = cell(nCond, nPart);  % RP per (cond, part)
RP_thresh_participants = cell(nCond, nPart);  % thresholded RP per (cond, part)
metrics_participants   = cell(nPart, 1);      % one metrics matrix [nCond x 8] per participant
metrics_tbl_participants = cell(nPart, 1);    % one table per participant

for part = 1:nPart
    sz = size(TimeSeries);
    non_singleton_dims = sum(sz > 1); %trick to get if the matrix is a vector
    if non_singleton_dims == 4
        disp(['Computing recurrence plots and RQA metrics for each participant.. ' num2str(part) ' / ' num2str(nPart)])
    else
        disp('Computing recurrence plots and RQA metrics..')
    end
    % --- RP per condition for this participant ---
    for cc = 1:nCond
        % Build RP for this (cond, participant)
        PS = phase_space_participants{cc, part};
        RP_temp_p        = zeros(size(PS,1));
        RP_thresh_temp_p = zeros(size(PS,1));
        for ii = 1:size(PS,1)        % over time-points
            for jj = 1:size(PS,1)    % over time-points
                RP_temp_p(ii,jj) = norm(PS(ii,:) - PS(jj,:));
            end
        end
        RP_thresh_temp_p(RP_temp_p < max(RP_temp_p(:)) * eps) = 1;
        if ~isempty(theiler_window)
            RP_thresh_temp_p(~theiler_mask) = 0; %exclude temporally adjacent recurrence points
        end

        RP_participants{cc, part}        = RP_temp_p;
        RP_thresh_participants{cc, part} = RP_thresh_temp_p;
        % clear RP_temp_p RP_thresh_temp_p
    end

    % --- RQA metrics for this participant (per condition) ---
    lmin = 2;
    vmin = 2;

    metrics_p = zeros(nCond, 8);
    for cc = 1:nCond
        RPth = RP_thresh_participants{cc, part};

        %%%%% RECURRENCE RATE (RR)
        if isempty(theiler_window)
            metrics_p(cc,1) = sum(RPth(:)) / numel(RPth);
        else
            metrics_p(cc,1) = sum(RPth(:)) / sum(theiler_mask(:));
        end

        %%%%% mean diagonal line (L) + collect diagonals
        [~, Ldiags] = dl(RPth);
        Ldiags(Ldiags < lmin) = [];
        metrics_p(cc,2) = mean(Ldiags);
        has_diagonal_lines = ~isempty(Ldiags);

        %%%%% Determinism (DET)
        if isempty(Ldiags)
            Ldiags = 0;
            warning('No diagonal elements in the RP (participant %d, condition %d)', part, cc);
        end
        if sum(RPth(:)) > 0
            metrics_p(cc,3) = sum(Ldiags) / sum(RPth(:));
        else
            metrics_p(cc,3) = NaN;
            warning('No recurrence points in the RP (participant %d, condition %d)', part, cc);
        end

        %%%%% Entropy (ENTR)
        histL = hist(Ldiags(:), 1:min(size(RPth)));
        metrics_p(cc,4) = entropy(histL(:));

        %%%%% Trapping time (TT) and Laminarity (LAM)
        [~, TTverts] = tt(RPth);
        TTverts(TTverts < vmin) = [];
        metrics_p(cc,5) = mean(TTverts);
        if sum(TTverts) > 0
            metrics_p(cc,6) = sum(TTverts) / sum(RPth(:));
        else
            metrics_p(cc,6) = NaN;
        end

        %%%%% Vmax
        if isempty(TTverts)
            metrics_p(cc,7) = NaN;
        else
            metrics_p(cc,7) = max(TTverts);
        end

        %%%%% DIV (inverse of maximum diagonal length)
        if isempty(theiler_window)
            if numel(Ldiags) >= 2
                Lmax = max(Ldiags(1:end-1));
                metrics_p(cc,8) = 1 / Lmax;
            elseif numel(Ldiags) == 1
                metrics_p(cc,8) = 1 / Ldiags(1);
            else
                metrics_p(cc,8) = NaN;
            end
        elseif has_diagonal_lines
            metrics_p(cc,8) = 1 / max(Ldiags);
        else
            metrics_p(cc,8) = NaN;
        end
    end

    headers = {'RR', 'L', 'DET', 'ENTR', 'TT', 'LAM', 'V_max', 'DIV'};
    metrics_tbl_participants{part} = array2table(metrics_p, 'VariableNames', headers);
    metrics_participants{part}     = metrics_p;
end

%% ---------------------------- Store outputs -----------------------------

% Phase-space coordinates are exposed so that they can be inspected or
% passed directly to BROADNESS_PhaseSpaceStatistics without recomputing
% the phase-space embedding or the recurrence analysis.
RQA_BROADNESS.PhaseSpace.Time                   = time(reduced_time_idx);
RQA_BROADNESS.PhaseSpace.PCs                    = PCs;
RQA_BROADNESS.PhaseSpace.ParticipantCoordinates = phase_space_participants;
RQA_BROADNESS.PhaseSpace.MeanCoordinates        = phase_space;
RQA_BROADNESS.PhaseSpace.nConditions            = nCond;
RQA_BROADNESS.PhaseSpace.nParticipants          = nPart;
RQA_BROADNESS.PhaseSpace.nDimensions            = length(PCs);
RQA_BROADNESS.PhaseSpace.Normalization.Method    = normalization;
RQA_BROADNESS.PhaseSpace.Normalization.Center    = normalizationCenter;
RQA_BROADNESS.PhaseSpace.Normalization.Scale     = normalizationScale;
RQA_BROADNESS.PhaseSpace.Normalization.Reference = ...
    'Selected time-points pooled across all conditions and participants';

% NEW: Per-participant recurrence plots + metrics
RQA_BROADNESS.RecurrencePlots.DistMat   = RP_participants;        % cell(nCond,nPart)
RQA_BROADNESS.RecurrencePlots.RecurPlot = RP_thresh_participants; % cell(nCond,nPart)
RQA_BROADNESS.RQA_metrics               = metrics_tbl_participants;  % {nPart} of tables
FIGURES.Handles = figureHandles;
FIGURES.Files = figureFiles;
RQA_BROADNESS.Figures.Files = figureFiles;

function plot_phase_space(ax, phaseSpace, analysedTime, condition)
nTimePoints = size(phaseSpace,1);
colors = jet(nTimePoints);
hold(ax, 'on'); grid(ax, 'minor'); box(ax, 'on');
if size(phaseSpace,2) == 2
    scatter(ax, phaseSpace(:,1), phaseSpace(:,2), 24, colors, 'filled');
    xlabel(ax, 'Brain network 1'); ylabel(ax, 'Brain network 2');
elseif size(phaseSpace,2) == 3
    scatter3(ax, phaseSpace(:,1), phaseSpace(:,2), phaseSpace(:,3), ...
        24, colors, 'filled');
    xlabel(ax, 'Brain network 1'); ylabel(ax, 'Brain network 2');
    zlabel(ax, 'Brain network 3'); view(ax, 3); axis(ax, 'vis3d');
else
    text(ax, 0.5, 0.5, 'Phase-space display requires 2 or 3 components', ...
        'HorizontalAlignment', 'center'); axis(ax, 'off');
    return
end
colormap(ax, jet);
c = colorbar(ax);
timeTicks = round(linspace(1, nTimePoints, min(6,nTimePoints)));
c.Ticks = linspace(0,1,numel(timeTicks));
c.TickLabels = round(analysedTime(timeTicks),3);
c.Label.String = 'Time (s)';
title(ax, ['Phase space — Condition ' num2str(condition)]);
axis(ax, 'square');
end

function plot_recurrence_matrix(ax, matrixToPlot, analysedTime, plotTitle, thresholded)
imagesc(ax, analysedTime, analysedTime, matrixToPlot);
set(ax, 'YDir', 'normal'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Time (s)');
c = colorbar(ax); axis(ax, 'square'); title(ax, plotTitle);
if thresholded
    colormap(ax, gray(2));
    caxis(ax, [0 1]);
    c.Ticks = [0 1];
    c.TickLabels = {'0','1'};
else
    colormap(ax, flipud(parula));
end
end

%% ------------------------ Helper: parse name/values ---------------------
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

end

