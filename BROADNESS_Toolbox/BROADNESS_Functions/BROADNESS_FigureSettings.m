function Settings = BROADNESS_FigureSettings(figureMode, figureLayout, outputPath, ...
    figureFormats, figurePrefix, functionFolder)

% ========================================================================
%  BROADNESS INTERNAL FIGURE SETTINGS
% ========================================================================
%
%  This internal utility validates and prepares the common figure options
%  used by the BROADNESS visualization functions.
%
%  figureMode   : 'off', 'show', 'save', or 'both'
%  figureLayout : 'individual', 'summary', or 'both'
%  outputPath   : Base output folder (required for 'save' and 'both')
%  figureFormats: Character vector, string, or cell array containing
%                 'png', 'pdf', and/or 'fig'
%  figurePrefix : Optional prefix added to every saved filename
%  functionFolder: Subfolder created inside OutputPath/BROADNESS_Figures
%
% ========================================================================

figureMode = lower(char(string(figureMode)));
figureLayout = lower(char(string(figureLayout)));

if ~ismember(figureMode, {'off','show','save','both'})
    error('"FigureMode" must be ''off'', ''show'', ''save'', or ''both''.')
end
if ~ismember(figureLayout, {'individual','summary','both'})
    error('"FigureLayout" must be ''individual'', ''summary'', or ''both''.')
end

if ischar(figureFormats) || (isstring(figureFormats) && isscalar(figureFormats))
    figureFormats = cellstr(figureFormats);
elseif isstring(figureFormats)
    figureFormats = cellstr(figureFormats(:));
elseif ~iscell(figureFormats)
    error('"FigureFormats" must be a character vector, string, or cell array.')
end

for formati = 1:length(figureFormats)
    figureFormats{formati} = lower(char(string(figureFormats{formati})));
end
figureFormats = unique(figureFormats, 'stable');
if isempty(figureFormats) || any(~ismember(figureFormats, {'png','pdf','fig'}))
    error('"FigureFormats" may contain only ''png'', ''pdf'', and ''fig''.')
end

if isempty(figurePrefix)
    figurePrefix = '';
else
    figurePrefix = matlab.lang.makeValidName(char(string(figurePrefix)));
    figurePrefix = [figurePrefix '_'];
end

Settings.Mode = figureMode;
Settings.Layout = figureLayout;
Settings.Formats = figureFormats;
Settings.Prefix = figurePrefix;
Settings.Show = ismember(figureMode, {'show','both'});
Settings.Save = ismember(figureMode, {'save','both'});
Settings.MakeIndividual = ismember(figureLayout, {'individual','both'});
Settings.MakeSummary = ismember(figureLayout, {'summary','both'});
Settings.Visible = 'off';
Settings.OpenFigureVisibility = 'invisible';
if Settings.Show
    Settings.Visible = 'on';
    Settings.OpenFigureVisibility = 'visible';
end
Settings.OutputFolder = '';

if Settings.Save
    if isempty(outputPath)
        error('An output path is required when "FigureMode" is ''save'' or ''both''.')
    end
    Settings.OutputFolder = fullfile(char(string(outputPath)), ...
        'BROADNESS_Figures', functionFolder);
    if ~exist(Settings.OutputFolder, 'dir')
        mkdir(Settings.OutputFolder);
    end
end

end
