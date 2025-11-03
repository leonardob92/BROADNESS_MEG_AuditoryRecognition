function [RQA_BROADNESS] = BROADNESS_PhaseSpace_RQA(BROADNESS, varargin)
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
%      - 'video'                        : 'on' or 'off' to show animated phase space plot (default: 'on')
%      - 'figure'                       : 'on' or 'off' to show the figures (default: 'on')
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - RQA_BROADNESS                              : Structure with RQA results
%      - .RecurrencePlots.DistMat               : Cell array of distance matrices
%      - .RecurrencePlots.RecurPlot             : Cell array of recurrence plots (i.e., thresholded distance matrices)
%      - .RQA_metrics                           : Table of 8 RQA measures
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
opts = struct('principalcomps', 1:2, 'timeinterval', [], 'threshold', 0.1, 'video', 'off','figure','off');
opts = parse_name_value_pairs(opts, varargin{:});

% Assign to readable internal names
PCs          = opts.principalcomps;
time_seconds = opts.timeinterval;
eps          = opts.threshold;
video        = opts.video;
figurel      = opts.figure;

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
% if ~(ismatrix(TimeSeries) || ndims(TimeSeries) == 3)
%     error(['"TimeSeries_BrainNetworks" must be 2D or 3D (time × components × [conditions]) ', ...
%            'or (components × time × [conditions]).']);
% end

% --- Validate optional args ---
if ~(isnumeric(PCs) && isvector(PCs))
    error('"principalcomps" must be a numeric vector.');
end
if ~isnumeric(eps) || ~isscalar(eps) || ~isfinite(eps)
    error('"threshold" must be a finite numeric scalar.');
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

% Indices for the chosen window (start >=, end <=)
timemin = find(time >= time_seconds(1), 1, 'first');
timemax = find(time <= time_seconds(2), 1, 'last');
reduced_time_idx = timemin:timemax;

isPCA = isfield(BROADNESS, 'Variance_BrainNetworks');  % PCA outputs variance field; ICA doesn't

if isPCA && length(reduced_time_idx) > size(TimeSeries,1)
    reduced_time_idx = reduced_time_idx(1:end-1); %adjusting potential mismatch of one timepoint between PCA time series and time..
    time = time(1:end-1);
end
    
if ~isPCA %reshaping matrix of time series if ICA
    TimeSeries = permute(TimeSeries,[2 1 3 4]); % 4th dimension if the time series were computed for each participant; this works even if the time series matrix is only 3D
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
   DM{cc} = RP_temp;
   RP{cc} = RP_thresh_temp;
%    clear RP_temp
%    clear RP_thresh_temp
end
    
%% ------------------------- Recurrence plots -----------------------------

if strcmp(figurel, 'on')
    % ------------------ Phase-space scatter (static) ------------------
    % Plot the averaged phase-space scatter plots (2D or 3D) so that
    % figure='on' produces the phase-space visuals even if video='off'.
    for cond = 1:numel(phase_space)
        PS = phase_space{cond};
        nTps = size(PS,1);
        cmap = jet(nTps);
        scatsize = 36;

        if size(PS,2) == 2
            figure;
            hold on;
            set(gcf,'Color','w');
            grid minor; box on;
            xlabel('Brain network 1'); ylabel('Brain network 2');
            title(['Phase space (static) - Condition ' num2str(cond)]);
            xlim([min(PS(:,1)) max(PS(:,1))]);
            ylim([min(PS(:,2)) max(PS(:,2))]);

            % color by time
            for ii = 1:nTps
                scatter(PS(ii,1), PS(ii,2), scatsize, cmap(ii,:), 'filled');
            end
            c = colorbar;
            colormap(jet);
            % create tick labels using the reduced_time_idx mapping to time
            t_ticks = round(linspace(1, nTps, min(10, nTps)));
            c.Ticks = linspace(0,1,numel(t_ticks));
            c.TickLabels = time(reduced_time_idx(t_ticks));
            hold off;
            axis square;

        elseif size(PS,2) == 3
            figure;
            hold on;
            set(gcf,'Color','w');
            view(3);
            axis vis3d;
            grid minor; box on;
            xlabel('Brain network 1'); ylabel('Brain network 2'); zlabel('Brain network 3');
            title(['Phase space (static, 3D) - Condition ' num2str(cond)]);
            xlim([min(PS(:,1)) max(PS(:,1))]);
            ylim([min(PS(:,2)) max(PS(:,2))]);
            zlim([min(PS(:,3)) max(PS(:,3))]);

            for ii = 1:nTps
                scatter3(PS(ii,1), PS(ii,2), PS(ii,3), scatsize, 'MarkerFaceColor', cmap(ii,:), 'MarkerEdgeColor', cmap(ii,:));
            end
            c = colorbar;
            colormap(jet);
            t_ticks = round(linspace(1, nTps, min(10, nTps)));
            c.Ticks = linspace(0,1,numel(t_ticks));
            c.TickLabels = time(reduced_time_idx(t_ticks));
            hold off;
            axis square;

        else
            % For >3D, skip static visualization (you already warn for video)
            warning('Static phase-space scatter is available only for 2 or 3 principal components.');
        end
    end

    % ----------------------- Recurrence plots -------------------------
    % Plotting the NO THRESHOLDED recurrence plot for each condition
    for cond = 1:size(DM,1) % over conditions
        figure;
        imagesc(time(reduced_time_idx), time(reduced_time_idx), DM{cond});
        xlabel('Time (s)'); ylabel('Time (s)'); set(gca,'YDir','normal');
        colorbar
        set(gcf,'Color','w')
        % Flip colormap: recurrence = blue, non-recurrence = yellow
        colormap(flipud(parula));
        title(['Condition ' num2str(cond) ' no thresholding'])
        axis square
    end

    % Plotting the THRESHOLDED recurrence plot for each condition
    for cond = 1:size(RP,1) % over conditions
        figure;
        imagesc(time(reduced_time_idx), time(reduced_time_idx), RP{cond});
        xlabel('Time (s)'); ylabel('Time (s)'); set(gca,'YDir','normal');
        colorbar
        set(gcf,'Color','w')
        title(['Condition ' num2str(cond) ' thresholding ' num2str(eps*100) '%'])
        axis square
    end
end



% if strcmp(figurel, 'on')
%     % Plotting the NO THRESHOLDED recurrence plot for each condition
%     for cond = 1:size(DM,1) % over conditions
%         figure;
%         imagesc(time(reduced_time_idx), time(reduced_time_idx), DM{cond}); 
%         xlabel('Time (s)'); ylabel('Time (s)'); set(gca,'YDir','normal');
%         colorbar
%         set(gcf,'Color','w')
%         % Flip colormap: recurrence = blue, non-recurrence = yellow
%         colormap(flipud(parula));  
%         title(['Condition ' num2str(cond) ' no thresholding'])
%         axis square
%     end
%     
%     % Plotting the THRESHOLDED recurrence plot for each condition
%     for cond = 1:size(RP,1) % over conditions
%         figure;
%         imagesc(time(reduced_time_idx), time(reduced_time_idx), RP{cond}); 
%         xlabel('Time (s)'); ylabel('Time (s)'); set(gca,'YDir','normal');
%         colorbar
%         set(gcf,'Color','w')
%         title(['Condition ' num2str(cond) ' thresholding ' num2str(eps*100) '%'])
%         axis square
%     end
% end

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
        metrics_p(cc,1) = sum(RPth(:)) / numel(RPth);

        %%%%% mean diagonal line (L) + collect diagonals
        [~, Ldiags] = dl(RPth);
        Ldiags(Ldiags < lmin) = [];
        metrics_p(cc,2) = mean(Ldiags);

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
        if numel(Ldiags) >= 2
            Lmax = max(Ldiags(1:end-1));
            metrics_p(cc,8) = 1 / Lmax;
        elseif numel(Ldiags) == 1
            metrics_p(cc,8) = 1 / Ldiags(1);
        else
            metrics_p(cc,8) = NaN;
        end
    end

    headers = {'RR', 'L', 'DET', 'ENTR', 'TT', 'LAM', 'V_max', 'DIV'};
    metrics_tbl_participants{part} = array2table(metrics_p, 'VariableNames', headers);
    metrics_participants{part}     = metrics_p;
end

%% ---------------------------- Store outputs -----------------------------

% NEW: Per-participant recurrence plots + metrics
RQA_BROADNESS.RecurrencePlots.DistMat   = RP_participants;        % cell(nCond,nPart)
RQA_BROADNESS.RecurrencePlots.RecurPlot = RP_thresh_participants; % cell(nCond,nPart)
RQA_BROADNESS.RQA_metrics               = metrics_tbl_participants;  % {nPart} of tables

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

