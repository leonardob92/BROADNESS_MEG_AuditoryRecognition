function BROADNESS_Plot_ActivationTimeseries(data, time, varargin)
%%
% =========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  TIME SERIES PLOTTING FUNCTION
% =========================================================================
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2025). Advanced Science. 
%  BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in Auditory Memory Networks.
%  https://doi.org/10.1002/advs.202507878
%
% =========================================================================
%
%
%  DESCRIPTION:
%  ------------
%  Plots the average time series of brain networks (or other signals) across
%  conditions and subjects, with optional standard error shading or dotted lines.
%
%  This function supports visualizing significant time windows using either:
%    - shaded background patches (default), or
%    - horizontal lines placed above or below the plotted data
%
%  Suitable for data in the BROADNESS toolbox format or any 4D matrix with:
%     [networks x time x subjects x conditions]
%
%  ------------------------------------------------------------------------
%  REQUIRED INPUTS:
%  ------------------------------------------------------------------------
%    data  : 4D matrix [networks x time x subjects x conditions]
%    time  : Vector of time points (in seconds)
%
%  ------------------------------------------------------------------------
%  OPTIONAL NAME-VALUE PAIRS:
%  ------------------------------------------------------------------------
%    'Groups'                   : Cell array of subject indices per group
%    'GroupLabels'              : Labels for subject groups (cell array)
%    'ConditionLabels'          : Cell array of condition names
%    'ConditionXGroupColors'    : Cell array of [nGroups x 3] RGB matrices (1 per condition)
%    'SignificantWindows'       : Nested cell array where each row defines one Y-level.
%                                 Each row is a cell array of [start end] time pairs.
%                                 Example:
%                                 {
%                                     {[0.4 0.6], [0.9 1.1]};       % Level 1
%                                     {[0.2 0.5]};                  % Level 2
%                                     {[0.3 0.4], [1.2 1.5]}        % Level 3
%                                 }
%    'SignificanceColors'       : [nLevels x 3] RGB matrix. One color per Y-level (row).
%                                 Each level uses a single color for all of its lines.
%    'SignificanceStyle'        : 'patch' (default) or 'line'
%    'SignificanceLinePosition' : 'above' (default) or 'below'
%    'SignificanceLineOffset'   : Distance from data to line (default = 0.01)
%    'SignificanceLinePadding'  : Space between line and plot frame (default = 0.05)
%    'YLimits'                  : Manual y-axis limits [min max]
%    'XLimits'                  : Manual x-axis limits [min max]
%    'STEStyle'                 : 1 = dotted SE, 2 = shaded SE (default = 2)
%    'Transparency'             : Alpha for shaded SE (default = 0.3)
%
%  ------------------------------------------------------------------------
%  OUTPUT:
%  ------------------------------------------------------------------------
%    One figure per network with condition averages and significance markings
%
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Leonardo Bonetti & Mattia Rosso
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  mattia.rosso@clin.au.dk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 22/07/2025
% ========================================================================









%% ----------------------------
% STEP 1: Set default options
% ----------------------------
opts = struct( ...
    'Groups', [], ...
    'GroupLabels', [], ...
    'ConditionLabels', [], ...
    'ConditionXGroupColors', [], ...
    'SignificantWindows', [], ...
    'SignificanceColors', [], ...
    'SignificanceStyle', 'patch', ...
    'SignificanceLinePosition', 'above', ...
    'SignificanceLineOffset', 0.01, ...
    'SignificanceLinePadding', 0.05, ...
    'YLimits', [], ...
    'XLimits', [], ...
    'STEStyle', 2, ...
    'Transparency', 0.3);

% Parse user-defined name-value pairs and override defaults
opts = parse_name_value_pairs(opts, varargin{:});


%% ----------------------------
% STEP 2: Validate inputs
% ----------------------------

% Ensure 'data' is a 4D numeric array
if ~isnumeric(data)
    error('Input "data" must be a numeric array: [networks x time x subjects x conditions].');
end
[nNet, nTime, nSubj, nConds] = size(data);  % Extract dimensions

% Warn if too many NaNs
nanRatio = sum(isnan(data(:))) / numel(data);
if nanRatio > 0.05
    warning('Data contains %.1f%% NaN values.', nanRatio * 100);
end

% Ensure 'time' is a numeric vector matching the time dimension
if ~isnumeric(time) || ~isvector(time)
    error('Input "time" must be a numeric vector.');
end
if length(time) ~= nTime
    error('Length of "time" must match the second dimension of "data".');
end

% Validate 'Groups' if provided
if ~isempty(opts.Groups)
    if ~iscell(opts.Groups)
        error('"Groups" must be a cell array of subject index vectors.');
    end
    for g = 1:length(opts.Groups)
        if any(opts.Groups{g} < 1 | opts.Groups{g} > nSubj)
            error('Invalid subject index in group %d.', g);
        end
        if isempty(opts.Groups{g})
            error('Group %d has no subjects.', g);
        end
    end
end

% Ensure group labels match group count
if ~isempty(opts.GroupLabels) && length(opts.GroupLabels) ~= length(opts.Groups)
    error('Length of "GroupLabels" must match number of groups.');
end

% Check if condition labels match number of conditions
if ~isempty(opts.ConditionLabels) && length(opts.ConditionLabels) ~= nConds
    warning('Length of "ConditionLabels" does not match number of conditions.');
end

% If no groups specified, assume one group of all subjects
if isempty(opts.Groups)
    opts.Groups = {1:nSubj};
end
nGroups = length(opts.Groups);

% Validate condition colors
if ~isempty(opts.ConditionXGroupColors)
    if ~iscell(opts.ConditionXGroupColors) || length(opts.ConditionXGroupColors) ~= nConds
        error('ConditionXGroupColors must be a cell array with one entry per condition.');
    end
    for c = 1:nConds
        colorMat = opts.ConditionXGroupColors{c};
        if size(colorMat,2) ~= 3 || size(colorMat,1) ~= length(opts.Groups)
            error('Color matrix for condition %d must be [%d x 3].', c, length(opts.Groups));
        end
        if any(colorMat(:) < 0 | colorMat(:) > 1)
            error('RGB values must be between 0 and 1.');
        end
    end
end

% Validate significant time windows and colors
if ~isempty(opts.SignificantWindows)
    if ~iscell(opts.SignificantWindows)
        error('"SignificantWindows" must be a cell array.');
    end
    
    % Validate significance windows with nested levels
    sigWinGroups = opts.SignificantWindows;
    if ~iscell(sigWinGroups) || ~iscell(sigWinGroups{1})
        error('"SignificantWindows" must be a nested cell array (cell array of cell arrays).');
    end
    
    for level = 1:length(sigWinGroups)
        levelWins = sigWinGroups{level};
        if ~iscell(levelWins)
            error('Each level in "SignificantWindows" must be a cell array of [start end] pairs.');
        end
        for w = 1:length(levelWins)
            win = levelWins{w};
            if ~isnumeric(win) || length(win) ~= 2 || win(1) > win(2)
                warning('Invalid window at level %d, index %d.', level, w);
            end
            if win(1) < min(time) || win(2) > max(time)
                warning('Window [%g %g] (level %d) is outside time range.', win(1), win(2), level);
            end
        end
    end
    
    if ~isempty(opts.SignificanceColors)
        if size(opts.SignificanceColors,1) < length(opts.SignificantWindows)
            warning('Fewer colors than significant windows.');
        end
        if any(opts.SignificanceColors(:) < 0 | opts.SignificanceColors(:) > 1)
            error('SignificanceColors must be RGB values between 0 and 1.');
        end
    end
end

% Ensure STE style is valid
if ~ismember(opts.STEStyle, [1 2])
    warning('"STEStyle" must be 1 or 2. Defaulting to 2.');
    opts.STEStyle = 2;
end

% Ensure transparency is valid
if ~isnumeric(opts.Transparency) || opts.Transparency < 0 || opts.Transparency > 1
    warning('Invalid "Transparency". Using default (0.3).');
    opts.Transparency = 0.3;
end

% Warn if STE will not be visible
if opts.STEStyle == 2 && opts.Transparency == 0
    warning('Shaded SE selected but transparency = 0.');
end

% Ensure valid significance line position
if ~ismember(lower(opts.SignificanceLinePosition), {'above','below'})
    warning('Invalid "SignificanceLinePosition". Using "above".');
    opts.SignificanceLinePosition = 'above';
end

% Warn if too many legend entries
if nConds * length(opts.Groups) > 20
    warning('Too many plotted lines (%d).', nConds * length(opts.Groups));
end


%% ----------------------------
% STEP 3: Prepare plotting specs
% ----------------------------

% Auto-generate missing condition labels
if isempty(opts.ConditionLabels)
    opts.ConditionLabels = {};
end
if length(opts.ConditionLabels) < nConds
    for i = (length(opts.ConditionLabels)+1):nConds
        opts.ConditionLabels{i} = ['Condition ' num2str(i)];
    end
end

% Auto-generate missing group labels
if isempty(opts.GroupLabels)
    for i = 1:nGroups
        opts.GroupLabels{i} = ['Group ' num2str(i)];
    end
end

% Default colors if none provided
if isempty(opts.ConditionXGroupColors)
    baseColors = lines(nConds);  % MATLAB line color palette
    opts.ConditionXGroupColors = cell(1, nConds);
    for c = 1:nConds
        col_main = baseColors(c,:);
        opts.ConditionXGroupColors{c} = [col_main; col_main*0.6];  % Default 2 groups
    end
end


%% ----------------------------
% STEP 4: Plotting
% ----------------------------
for net = 1:nNet
    figure; hold on;
%     title(['Brain Network ' num2str(net)]);  % Title for each network
    
    all_vals = [];  % Store all y-values for axis scaling
    
    % Loop through conditions and groups
    for c = 1:nConds
        for g = 1:nGroups
            subj_idx = opts.Groups{g};  % Subjects in this group
            netdata = squeeze(data(net,:,subj_idx,c));  % time x subjects matrix
            mean_ts = mean(netdata, 2, 'omitnan');       % Mean across subjects
            ste_ts = std(netdata, 0, 2, 'omitnan') ./ sqrt(size(netdata,2));  % Standard Error
            
            color = opts.ConditionXGroupColors{c}(g,:);  % Color for this line
            
            % Compose legend label based on number of groups
            if nGroups == 1
                legendLabel = opts.ConditionLabels{c};
            else
                legendLabel = [opts.ConditionLabels{c} ' - ' opts.GroupLabels{g}];
            end
            
            % Plot mean line
            plot(time, mean_ts, 'Color', color, 'LineWidth', 1.5, ...
                'DisplayName', legendLabel);
            
            % Add standard error (STE)
            if opts.STEStyle == 1  % Dotted lines
                plot(time, mean_ts + ste_ts, ':', 'Color', color, 'HandleVisibility','off');
                plot(time, mean_ts - ste_ts, ':', 'Color', color, 'HandleVisibility','off');
            elseif opts.STEStyle == 2  % Shaded area
                x_fill = [time(:)' fliplr(time(:)')];
                y_fill = [(mean_ts + ste_ts)' fliplr((mean_ts - ste_ts)')];
                fill(x_fill, y_fill, color, 'FaceAlpha', opts.Transparency, 'EdgeColor','none', 'HandleVisibility','off');
            end
            
            % Store values for axis scaling
            all_vals = [all_vals; mean_ts + ste_ts; mean_ts - ste_ts];
        end
    end
    
    % Set Y-limits
    if isempty(opts.YLimits)
        margin = 0.05 * range(all_vals);
        ylims = [min(all_vals)-margin, max(all_vals)+margin];
        ylim(ylims);
    else
        ylims = opts.YLimits;
        ylim(ylims);
    end
    
    % Add significance markings (grouped by Y-levels, same color per level)
    if ~isempty(opts.SignificantWindows)
        sigWinGroups = opts.SignificantWindows;
        
        if ~iscell(sigWinGroups) || ~iscell(sigWinGroups{1})
            error('"SignificantWindows" must be a nested cell array: each row = one Y-level.');
        end
        
        nLevels = length(sigWinGroups);
        
        y_range = diff(ylims);
        offset = opts.SignificanceLineOffset * y_range;
        padding = opts.SignificanceLinePadding * y_range;
        line_spacing = offset;
        
        % Total height to fit all stacked levels
        total_height = (nLevels - 1) * line_spacing + padding;
        
        switch lower(opts.SignificanceLinePosition)
            case 'above'
                base_y = ylims(2) + offset;
                new_ylim = [ylims(1), base_y + total_height];
            case 'below'
                base_y = ylims(1) - offset;
                new_ylim = [base_y - total_height, ylims(2)];
            otherwise
                warning('Invalid SignificanceLinePosition. Defaulting to "above".');
                base_y = ylims(2) + offset;
                new_ylim = [ylims(1), base_y + total_height];
        end
        
        ylim(new_ylim);  % Expand Y-limits
        
        % Plot each level of significance windows
        for level = 1:nLevels
            levelWins = sigWinGroups{level};
            if ~iscell(levelWins)
                error('Each level in "SignificantWindows" must be a cell array of [start end] pairs.');
            end
            
            % Determine Y-coordinate for this level
            switch lower(opts.SignificanceLinePosition)
                case 'above'
                    y_sig = base_y + (level - 1) * line_spacing;
                case 'below'
                    y_sig = base_y - (level - 1) * line_spacing;
            end
            
            % Use same color for all lines at this level
            if isempty(opts.SignificanceColors) || size(opts.SignificanceColors,1) < level
                col = [0.85 0.85 0.85];  % Default gray
            else
                col = opts.SignificanceColors(level,:);
            end
            
            for w = 1:length(levelWins)
                win = levelWins{w};
                
                switch lower(opts.SignificanceStyle)
                    case 'patch'
                        patch([win(1) win(2) win(2) win(1)], [-1e5 -1e5 1e5 1e5], col, ...
                            'EdgeColor','none','FaceAlpha',0.2, 'HandleVisibility','off');
                    case 'line'
                        line([win(1) win(2)], [y_sig y_sig], 'Color', col, 'LineWidth', 2, 'HandleVisibility','off');
                    otherwise
                        warning('Unknown SignificanceStyle. Using patch.');
                        patch([win(1) win(2) win(2) win(1)], [-1e5 -1e5 1e5 1e5], col, ...
                            'EdgeColor','none','FaceAlpha',0.2, 'HandleVisibility','off');
                end
            end
        end
    end
    
    % Set X-limits
    if isempty(opts.XLimits)
        xlim([min(time), max(time)]);
    else
        xlim(opts.XLimits);
    end
    
    % Label axes and finalize plot
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Location','best');
    set(gcf,'Color','w');
    box on;
    grid minor;
end
end

%% ----------------------------
% Helper Function: parse_name_value_pairs
% ----------------------------
function opts = parse_name_value_pairs(opts, varargin)
% Parses name-value pair arguments and updates 'opts' structure
if mod(length(varargin), 2) ~= 0
    error('Arguments must be name-value pairs.');
end
validFields = fieldnames(opts);
for i = 1:2:length(varargin)
    name = varargin{i};
    value = varargin{i+1};
    match = strcmpi(name, validFields);  % Case-insensitive match
    if any(match)
        opts.(validFields{match}) = value;  % Update option
    else
        error(['Unknown parameter: ' name]);
    end
end
end


