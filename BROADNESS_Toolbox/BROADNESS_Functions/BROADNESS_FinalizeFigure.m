function savedFiles = BROADNESS_FinalizeFigure(fig, Settings, fileName)

% ========================================================================
%  BROADNESS INTERNAL FIGURE EXPORT
% ========================================================================
%
%  Applies common publication-oriented formatting, exports the requested
%  file formats, and closes figures generated in save-only mode.
%
% ========================================================================

savedFiles = {};
if isempty(fig) || ~isgraphics(fig, 'figure')
    return
end

set(fig, 'Color', 'w');
fontObjects = findall(fig, '-property', 'FontName');
if ~isempty(fontObjects)
    set(fontObjects, 'FontName', 'Helvetica Neue');
end

if Settings.Save
    fileName = matlab.lang.makeValidName(char(string(fileName)));
    for formati = 1:length(Settings.Formats)
        fileFormat = Settings.Formats{formati};
        outputFile = fullfile(Settings.OutputFolder, ...
            [Settings.Prefix fileName '.' fileFormat]);
        switch fileFormat
            case 'png'
                exportgraphics(fig, outputFile, 'Resolution', 300);
            case 'pdf'
                exportgraphics(fig, outputFile, 'ContentType', 'vector');
            case 'fig'
                savefig(fig, outputFile, 'compact');
        end
        savedFiles{end+1,1} = outputFile; %#ok<AGROW>
    end
end

if ~Settings.Show
    close(fig)
end

end
