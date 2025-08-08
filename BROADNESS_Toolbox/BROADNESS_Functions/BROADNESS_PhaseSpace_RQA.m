function [RQA_BROADNESS] = BROADNESS_PhaseSpace_RQA(BROADNESS, varargin)
%%
% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  RECURRENCE QUANTIFICATION ANALYSIS (RQA)
% ========================================================================
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P., Kringelbach, M.L., & Rosso, M. (2024).
%  BROADband brain Network Estimation via Source Separation (BROAD-NESS). bioRxiv, 2024-10.
%  https://doi.org/10.1101/2024.10.31.621257
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
%      - .TimeSeries_BrainNetworks      : 2D or 3D matrix (time × components × [conditions])
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
%      - .RQA_matrices.RecurrencePlot_NOthresh  : Cell array of recurrence plots (distance matrices)
%      - .RQA_matrices.RecurrencePlot_thresh    : Cell array of thresholded recurrence plots
%      - .RQA_metrics                           : Table of 8 RQA measures per experimental condition
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
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 08/08/2025
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
opts = struct('principalcomps', 1:10, 'timeinterval', [], 'threshold', 0.1, 'video', 'on','figure','on');
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
if ~(ismatrix(TimeSeries) || ndims(TimeSeries) == 3)
    error(['"TimeSeries_BrainNetworks" must be 2D or 3D (time × components × [conditions]) ', ...
           'or (components × time × [conditions]).']);
end

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
    TimeSeries = permute(TimeSeries,[2 1 3]);
end

%% ------------------- Compute phase space coordinates --------------------

% initializating the cell array that will contain the coordinates in the
% phase space for each condition
phase_space = cell(size(TimeSeries,3),1);
for cond = 1:size(TimeSeries,3) % over conditions

    % initializating the matrix that will contain the coordinates in the phase space
    phase_space_temp = zeros(length(reduced_time_idx),length(PCs));

        for t = 1:length(reduced_time_idx) % over time points
            for cc = 1:length(PCs) % over components
                phase_space_temp(t,cc) = TimeSeries(reduced_time_idx(t),PCs(cc),cond);
            end
        end
    phase_space{cond} = phase_space_temp;
%     clear phase_space_temp
        
end
    
%% ---------------------- Display phase space video -----------------------

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

% Initializing the cell array that will contain the recurrence plot for each condition 
RP = cell(size(phase_space,1), 1);

% Initializing the cell array that will contain the thresholded recurrence plot for each condition 
RP_thresh = cell(size(phase_space,1), 1);


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
   RP{cc} = RP_temp;
   RP_thresh{cc} =RP_thresh_temp;
%    clear RP_temp
%    clear RP_thresh_temp
   
end
    
%% ------------------------- Recurrence plots -----------------------------

if strcmp(figurel, 'on')
    % Plotting the NO THRESHOLDED recurrence plot for each condition
    for cond = 1:size(RP,1) % over conditions
        figure;
        imagesc(time(reduced_time_idx),time(reduced_time_idx), RP{cond}); xlabel('Time (s)'); ylabel('Time (s)'); set(gca,'YDir','normal');
        colorbar
        set(gcf,'Color','w')
        % caxis([0.8 2.8])
        title(['Condition ' num2str(cond) ' no thresholding'])
    end
    
    % Plotting the THRESHOLDED recurrence plot for each condition
    for cond = 1:size(RP_thresh,1) % over conditions
        figure;
        imagesc(time(reduced_time_idx),time(reduced_time_idx), RP_thresh{cond}); xlabel('Time (s)'); ylabel('Time (s)'); set(gca,'YDir','normal');
        colorbar
        set(gcf,'Color','w')
        title(['Condition ' num2str(cond) ' thresholding' num2str(eps*100) '%'])
    end
end

%% --------------------------- Computing metrics --------------------------

lmin = 2; % minimum consecutive elements on diagonal lines

vmin = 2; % minimum consecutive elements on vertical lines

% Initializing the matrix that will contain all the metrics
metrics = zeros(size(RP_thresh,1), 8);

for cc = 1:size(RP_thresh,1) % over conditions
    
    disp(['Computing 8 metrics - condition ' num2str(cc)])
    %%%%%%%%%%% RECURRENCE RATE (RR) %%%%%%%%%
    %number of recurrence point divided by total number of points
    metrics(cc,1) = sum(RP_thresh{cc}(:)) / numel(RP_thresh{cc});

    %%%%%%%%%%% mean diagonal line (L) %%%%%%%%%%%%%%%%
%     clear Ldiags
    %Ldiags contains length of the diagonals found in the plot (showing that something is recurrent (with a time-lag)
    [~ , Ldiags] = dl(RP_thresh{cc});
    Ldiags(Ldiags<lmin) = [];
    %Storing the mean of diagonals length
    metrics(cc, 2) = mean(Ldiags);

    %%%%%%%%%%% Determinism (DET) %%%%%%%%%%%%%%%
    if isempty(Ldiags)
        Ldiags = 0;
        warning('No diagonal elements in the RP');
    end

    if sum(RP_thresh{cc}(:)) > 0
        metrics(cc,3) = sum(Ldiags) / (sum(RP_thresh{cc}(:))); %all diagonal summed divided by total number of recurrences = proportion of elements arranged in diagonals compared to all elements
    else
        metrics(cc,3) = NaN;
        warning('No recurrence points in the RP');
    end
    %%%%%%%%%%%% Entropy (ENTR) %%%%%%%%%%%%%%%
     %distribution of the number of occurrences of the diagonal of different lengths
    histL = hist(Ldiags(:), 1:min(size(RP_thresh{cc})));
    metrics(cc,4) = entropy(histL(:));

    %%%%%%%%%%%% Trapping time (TT) %%%%%%%%%%%%%
%     clear TTverts
    %TTverts contains legnth of the vertical lines in the plot (showing that the plot is trapped in a state)
    [~, TTverts] = tt(RP_thresh{cc});
    TTverts(TTverts < vmin) = [];
    %Storing the mean vertical lines length
    metrics(cc,5) = mean(TTverts);

    %%%%%%%%%%%% Laminarity (LAM) %%%%%%%%%%%%%%%
    if sum(TTverts)>0
        %all vertical lines summed divided by total number of recurrences = proportion of elements arranged in vertical lines (being trapped in the same state) compared to all elements
        metrics(cc,6) = sum(TTverts) / (sum(RP_thresh{cc}(:)));
    else
      metrics(cc,6) = NaN;
    end
     
    %%%%%%%%%%% Maximal duration of laminar state, i.e, vertical line (Vmax) %%%%%%%%%%% 
    
    metrics(cc,7) = max(TTverts);
    %%%%%%%%%%% Divergence, i.e. inverse of maximum diagonal length (DIV) %%%%%%%%%%%
%     clear Lmax
    Lmax = max(Ldiags(1:end-1));
    metrics(cc,8) = 1 / Lmax;
    
end
headers = {'RR', 'L', 'DET', 'ENTR', 'TT', 'LAM', 'V_max', 'DIV'};
metrics_tbl = array2table(metrics, 'VariableNames', headers);

%% ---------------------------- Store outputs -----------------------------

RQA_BROADNESS.RecurrencePlots.RecurrencePlot_NOthresh = RP;
RQA_BROADNESS.RecurrencePlots.RecurrencePlot_thresh = RP_thresh;
RQA_BROADNESS.RQA_metrics = metrics_tbl;

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