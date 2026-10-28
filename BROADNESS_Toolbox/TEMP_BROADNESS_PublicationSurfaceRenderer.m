function [fig, RENDER] = TEMP_BROADNESS_PublicationSurfaceRenderer( ...
    spatialPattern, MNI_coords, varargin)
%%
% ========================================================================
%  BROADNESS EXPERIMENTAL PUBLICATION SURFACE RENDERER
% ========================================================================
%
%  This experimental function projects one source-space spatial activation
%  pattern onto the full-brain MNI152 surface bundled with BROADNESS. Unlike
%  the standard dot visualization, it generates a continuous activation
%  overlay on the cortical and cerebellar surface.
%
%  Example:
%       load('MNI152_8mm_coord_dyi.mat','MNI8')
%       spatialPattern = BROADNESS.ActivationPatterns_BrainNetworks(:,1);
%       TEMP_BROADNESS_PublicationSurfaceRenderer(spatialPattern,MNI8);
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - spatialPattern       : One value for every brain source (sources x 1)
%  - MNI_coords           : MNI coordinates in the same order (sources x 3)
%
%  Optional name-value arguments:
%  - 'threshold'          : Absolute surface threshold. By default it is
%                           mean(abs(values)) + thresholdsd*std(abs(values))
%  - 'thresholdsd'        : Standard deviations used by the automatic
%                           threshold (default: 1)
%  - 'projectionradius'   : Maximum source-to-surface projection distance
%                           in millimetres (default: 18)
%  - 'kernelwidth'        : Gaussian projection width in mm (default: 6)
%  - 'smoothingiterations': Surface-smoothing iterations (default: 2)
%  - 'smoothingstrength'  : Neighbour contribution per iteration, between
%                           0 and 1 (default: 0.25)
%  - 'title'              : Figure title (default: 'Spatial activation pattern')
%  - 'visible'            : 'on' or 'off' (default: 'on')
%  - 'outputfile'         : Optional .png, .pdf, or .fig output filename
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - fig                  : MATLAB figure handle
%  - RENDER               : Projected surface values, threshold, surface
%                           geometry, and rendering settings
%
% ------------------------------------------------------------------------
%  IMPORTANT:
% ------------------------------------------------------------------------
%  - This function improves visualization only. It does not increase the
%    spatial resolution of the original source data.
%  - Statistical testing must be performed on the original source values,
%    not on the interpolated surface values produced here.
%  - Positive and negative weights retain their original polarity. Their
%    interpretation follows the orientation chosen for the PCA component.
%
% ========================================================================

%% ----------------------------- Inputs ----------------------------------

opts = struct( ...
    'threshold', [], ...
    'thresholdsd', 1, ...
    'projectionradius', 18, ...
    'kernelwidth', 6, ...
    'smoothingiterations', 2, ...
    'smoothingstrength', 0.25, ...
    'title', 'Spatial activation pattern', ...
    'visible', 'on', ...
    'outputfile', '');
opts = parse_name_value_pairs(opts,varargin{:});

spatialPattern = spatialPattern(:);
if ~isnumeric(spatialPattern) || ~isreal(spatialPattern) || isempty(spatialPattern)
    error('"spatialPattern" must be a non-empty real numeric vector.');
end
if ~isnumeric(MNI_coords) || size(MNI_coords,2) ~= 3 || ...
        size(MNI_coords,1) ~= length(spatialPattern)
    error('"MNI_coords" must contain one row per spatial-pattern value.');
end
validSources = isfinite(spatialPattern) & all(isfinite(MNI_coords),2);
if sum(validSources) < 4
    error('At least four sources must contain finite values and coordinates.');
end

validate_positive_scalar(opts.thresholdsd,'thresholdsd',true)
validate_positive_scalar(opts.projectionradius,'projectionradius',false)
validate_positive_scalar(opts.kernelwidth,'kernelwidth',false)
if ~isnumeric(opts.smoothingiterations) || ~isscalar(opts.smoothingiterations) || ...
        ~isfinite(opts.smoothingiterations) || opts.smoothingiterations < 0 || ...
        fix(opts.smoothingiterations) ~= opts.smoothingiterations
    error('"smoothingiterations" must be a non-negative integer.');
end
if ~isnumeric(opts.smoothingstrength) || ~isscalar(opts.smoothingstrength) || ...
        ~isfinite(opts.smoothingstrength) || opts.smoothingstrength < 0 || ...
        opts.smoothingstrength > 1
    error('"smoothingstrength" must be between 0 and 1.');
end
if ~isempty(opts.threshold) && (~isnumeric(opts.threshold) || ...
        ~isscalar(opts.threshold) || ~isfinite(opts.threshold) || ...
        opts.threshold < 0)
    error('"threshold" must be empty or a non-negative numeric scalar.');
end
if ~ismember(lower(char(string(opts.visible))),{'on','off'})
    error('"visible" must be ''on'' or ''off''.');
end

templateFile = which('BrainTemplate_MNI152_1mm_FullBrain.fig');
if isempty(templateFile)
    localTemplate = fullfile(fileparts(mfilename('fullpath')), ...
        'BROADNESS_External','BrainTemplate_MNI152_1mm_FullBrain.fig');
    if exist(localTemplate,'file')
        templateFile = localTemplate;
    else
        error(['The full-brain template was not found. Run BROADNESS_Startup ' ...
            'or restore BrainTemplate_MNI152_1mm_FullBrain.fig.']);
    end
end

%% --------------------- Obtain the surface geometry ---------------------

templateFigure = openfig(templateFile,'new','invisible');
templatePatch = findobj(templateFigure,'Type','patch');
if isempty(templatePatch)
    close(templateFigure)
    error('The full-brain template does not contain a surface patch.');
end
vertices = templatePatch(1).Vertices;
faces = templatePatch(1).Faces;
close(templateFigure)

%% ------------------ Project source values to the surface ----------------

sourceCoordinates = MNI_coords(validSources,:);
sourceValues = spatialPattern(validSources);
surfaceValues = gaussian_surface_projection(vertices,sourceCoordinates, ...
    sourceValues,opts.projectionradius,opts.kernelwidth);

if opts.smoothingiterations > 0
    surfaceValues = smooth_surface_values(surfaceValues,faces, ...
        opts.smoothingiterations,opts.smoothingstrength);
end

finiteSurface = isfinite(surfaceValues);
if ~any(finiteSurface)
    error(['No source was sufficiently close to the template surface. ' ...
        'Increase "projectionradius".']);
end

if isempty(opts.threshold)
    absoluteValues = abs(surfaceValues(finiteSurface));
    activationThreshold = mean(absoluteValues) + ...
        opts.thresholdsd*std(absoluteValues);
else
    activationThreshold = opts.threshold;
end
maximumAbsoluteValue = max(abs(surfaceValues(finiteSurface)));
if maximumAbsoluteValue == 0
    maximumAbsoluteValue = 1;
end

% A high, nearly constant colour contribution above threshold creates
% clearly bounded activation areas without changing the underlying values.
surfaceAlpha = zeros(size(surfaceValues));
aboveThreshold = finiteSurface & abs(surfaceValues) >= activationThreshold;
if maximumAbsoluteValue > activationThreshold
    surfaceAlpha(aboveThreshold) = 0.88 + 0.12 * ...
        (abs(surfaceValues(aboveThreshold))-activationThreshold) ./ ...
        (maximumAbsoluteValue-activationThreshold);
else
    surfaceAlpha(aboveThreshold) = 1;
end
surfaceAlpha = min(1,max(0,surfaceAlpha));
surfaceColorMap = publication_colormap(256);
surfaceColors = blend_surface_colors(surfaceValues,surfaceAlpha, ...
    maximumAbsoluteValue,surfaceColorMap);

%% -------------------------- Four-view figure ----------------------------

fig = figure('Visible',lower(char(string(opts.visible))),'Color','w', ...
    'Position',[70 70 1100 760],'Name',char(string(opts.title)), ...
    'NumberTitle','off');
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

faceCentresX = mean(reshape(vertices(faces,1),size(faces)),2);
leftFaces = faces(faceCentresX <= 1,:);
rightFaces = faces(faceCentresX >= -1,:);

viewDefinitions = { ...
    leftFaces,  [-90 8], 'Left lateral'; ...
    rightFaces, [90 8],  'Right lateral'; ...
    leftFaces,  [90 8],  'Left medial'; ...
    rightFaces, [-90 8], 'Right medial'};

brainAxes = gobjects(size(viewDefinitions,1),1);
for viewIndex = 1:size(viewDefinitions,1)
    ax = nexttile(layout);
    brainAxes(viewIndex) = ax;
    draw_surface_view(ax,vertices,viewDefinitions{viewIndex,1}, ...
        surfaceColors,viewDefinitions{viewIndex,2});
end

colormap(fig,surfaceColorMap)
set(findall(fig,'Type','axes'),'CLim',[-maximumAbsoluteValue maximumAbsoluteValue])
title(layout,char(string(opts.title)),'FontWeight','normal','FontSize',17)

% A dedicated invisible axis provides one shared colour bar for all views.
colorAxis = axes(fig,'Position',[0.93 0.18 0.015 0.64],'Visible','off', ...
    'CLim',[-maximumAbsoluteValue maximumAbsoluteValue]);
colormap(colorAxis,surfaceColorMap)
colorBar = colorbar(colorAxis,'Location','eastoutside');
colorBar.Label.String = 'Spatial activation weight';
colorBar.FontName = 'Helvetica Neue';
colorBar.FontSize = 11;

% Figure annotations keep view labels in a stable position for both screen
% rendering and exported graphics, including perspective 3D axes.
drawnow
for viewIndex = 1:length(brainAxes)
    axesPosition = brainAxes(viewIndex).Position;
    annotation(fig,'textbox',[axesPosition(1), ...
        axesPosition(2)+axesPosition(4)-0.035,axesPosition(3),0.03], ...
        'String',viewDefinitions{viewIndex,3},'LineStyle','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Helvetica Neue','FontSize',13);
end

fontObjects = findall(fig,'-property','FontName');
if ~isempty(fontObjects)
    set(fontObjects,'FontName','Helvetica Neue')
end

%% ----------------------------- Output ----------------------------------

RENDER.Vertices = vertices;
RENDER.Faces = faces;
RENDER.SurfaceValues = surfaceValues;
RENDER.SurfaceAlpha = surfaceAlpha;
RENDER.Threshold = activationThreshold;
RENDER.MaximumAbsoluteValue = maximumAbsoluteValue;
RENDER.ValidSources = validSources;
RENDER.Settings = opts;
RENDER.Interpretation = ['Surface values are Gaussian projections for ' ...
    'visualization only; inference remains in the original source space.'];

if ~isempty(opts.outputfile)
    outputFile = char(string(opts.outputfile));
    [outputFolder,~,extension] = fileparts(outputFile);
    if ~isempty(outputFolder) && ~exist(outputFolder,'dir')
        mkdir(outputFolder)
    end
    switch lower(extension)
        case '.png'
            exportgraphics(fig,outputFile,'Resolution',300)
        case '.pdf'
            exportgraphics(fig,outputFile,'ContentType','vector')
        case '.fig'
            savefig(fig,outputFile,'compact')
        otherwise
            error('"outputfile" must end in .png, .pdf, or .fig.');
    end
    RENDER.OutputFile = outputFile;
else
    RENDER.OutputFile = '';
end

end


%% ========================================================================
%                              Helper functions
% =========================================================================

function surfaceValues = gaussian_surface_projection(vertices,coordinates, ...
        values,projectionRadius,kernelWidth)
numberVertices = size(vertices,1);
surfaceValues = nan(numberVertices,1);
coordinateNorms = sum(coordinates.^2,2)';
radiusSquared = projectionRadius^2;
blockSize = 250;

for firstVertex = 1:blockSize:numberVertices
    lastVertex = min(numberVertices,firstVertex+blockSize-1);
    blockVertices = vertices(firstVertex:lastVertex,:);
    squaredDistances = bsxfun(@plus,sum(blockVertices.^2,2),coordinateNorms) - ...
        2*(blockVertices*coordinates');
    squaredDistances = max(squaredDistances,0);
    withinRadius = squaredDistances <= radiusSquared;
    weights = exp(-squaredDistances/(2*kernelWidth^2)) .* withinRadius;
    weightSums = sum(weights,2);
    projected = (weights*values) ./ weightSums;
    projected(weightSums == 0) = NaN;
    surfaceValues(firstVertex:lastVertex) = projected;
end
end


function values = smooth_surface_values(values,faces,iterations,strength)
numberVertices = length(values);
rows = [faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3); faces(:,1)];
columns = [faces(:,2); faces(:,1); faces(:,3); faces(:,2); faces(:,1); faces(:,3)];
adjacency = sparse(rows,columns,1,numberVertices,numberVertices);
adjacency = adjacency > 0;

for iteration = 1:iterations
    finiteValues = isfinite(values);
    valuesForAverage = values;
    valuesForAverage(~finiteValues) = 0;
    weightedValues = adjacency*valuesForAverage;
    neighbourCounts = adjacency*double(finiteValues);
    neighbourMean = weightedValues ./ neighbourCounts;
    update = neighbourCounts > 0;
    values(update) = (1-strength)*values(update) + ...
        strength*neighbourMean(update);
end
end


function draw_surface_view(ax,vertices,faces,vertexColors,cameraView)
hold(ax,'on')

% A single opaque surface blends the anatomical grey with activation
% colours. This avoids transparency-sorting artefacts during PDF/PNG export.
patch(ax,'Vertices',vertices,'Faces',faces, ...
    'FaceVertexCData',vertexColors,'FaceColor','interp','FaceAlpha',1, ...
    'EdgeColor','none','FaceLighting','gouraud', ...
    'AmbientStrength',0.48,'DiffuseStrength',0.72, ...
    'SpecularStrength',0.04,'BackFaceLighting','reverselit');

plot_mesh_silhouette(ax,vertices,faces,cameraView)

axis(ax,'equal')
axis(ax,'vis3d')
axis(ax,'off')
view(ax,cameraView)
camzoom(ax,1.75)
camlight(ax,'headlight')
lighting(ax,'gouraud')
material(ax,'dull')
end


function plot_mesh_silhouette(ax,vertices,faces,cameraView)
% Identify edges separating camera-facing and rear-facing triangles. These
% lines provide a clean anatomical outline without displaying the full mesh.
firstVector = vertices(faces(:,2),:)-vertices(faces(:,1),:);
secondVector = vertices(faces(:,3),:)-vertices(faces(:,1),:);
faceNormals = cross(firstVector,secondVector,2);
azimuth = cameraView(1);
elevation = cameraView(2);
cameraDirection = [sind(azimuth)*cosd(elevation), ...
    -cosd(azimuth)*cosd(elevation),sind(elevation)];
cameraFacing = faceNormals*cameraDirection' >= 0;

allEdges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
allEdges = sort(allEdges,2);
faceIndices = repmat((1:size(faces,1))',3,1);
[uniqueEdges,~,edgeGroups] = unique(allEdges,'rows');
edgeCounts = accumarray(edgeGroups,1);
facingCounts = accumarray(edgeGroups,double(cameraFacing(faceIndices)));
outlineMask = edgeCounts == 1 | ...
    (facingCounts > 0 & facingCounts < edgeCounts);
outlineEdges = uniqueEdges(outlineMask,:);

if isempty(outlineEdges)
    return
end

% Remove very short disconnected fragments, which otherwise look like
% isolated strokes rather than meaningful anatomical contours.
outlineGraph = graph(outlineEdges(:,1),outlineEdges(:,2),[],size(vertices,1));
outlineComponents = conncomp(outlineGraph);
edgeComponents = outlineComponents(outlineEdges(:,1))';
componentSizes = accumarray(edgeComponents,1);
outlineEdges = outlineEdges(componentSizes(edgeComponents) >= 5,:);
if isempty(outlineEdges)
    return
end

centre = mean(vertices,1);
lineVertices = bsxfun(@plus,centre,1.002*bsxfun(@minus,vertices,centre));
x = [lineVertices(outlineEdges(:,1),1)'; ...
    lineVertices(outlineEdges(:,2),1)'; nan(1,size(outlineEdges,1))];
y = [lineVertices(outlineEdges(:,1),2)'; ...
    lineVertices(outlineEdges(:,2),2)'; nan(1,size(outlineEdges,1))];
z = [lineVertices(outlineEdges(:,1),3)'; ...
    lineVertices(outlineEdges(:,2),3)'; nan(1,size(outlineEdges,1))];
plot3(ax,x(:),y(:),z(:),'-','Color',[0.22 0.23 0.25], ...
    'LineWidth',0.55,'HandleVisibility','off')
end


function colors = blend_surface_colors(values,alphaValues,maximumValue,colorMap)
baseColor = [0.60 0.62 0.65];
numberColors = size(colorMap,1);
normalizedValues = (values+maximumValue)/(2*maximumValue);
normalizedValues = min(1,max(0,normalizedValues));
colorIndices = round(normalizedValues*(numberColors-1))+1;
colorIndices(~isfinite(colorIndices)) = ceil(numberColors/2);
activationColors = colorMap(colorIndices,:);
colors = bsxfun(@times,1-alphaValues,baseColor) + ...
    bsxfun(@times,alphaValues,activationColors);
colors(~isfinite(values),:) = repmat(baseColor,sum(~isfinite(values)),1);
end


function colors = publication_colormap(numberColors)
numberNegative = floor(numberColors/2);
numberPositive = numberColors-numberNegative;

negativeStart = [0.05 0.13 0.72];
negativeEnd = [0.10 0.32 0.92];
positiveStart = [0.68 0.02 0.01];
positiveEnd = [1.00 0.88 0.02];

negativeColors = interpolate_colors(negativeStart,negativeEnd,numberNegative);
positiveColors = interpolate_colors(positiveStart,positiveEnd,numberPositive);
colors = [negativeColors; positiveColors];
end


function colors = interpolate_colors(firstColor,lastColor,numberColors)
mixing = linspace(0,1,numberColors)';
colors = bsxfun(@times,1-mixing,firstColor) + ...
    bsxfun(@times,mixing,lastColor);
end


function validate_positive_scalar(value,name,allowZero)
lowerLimit = 0;
if allowZero
    invalid = value < lowerLimit;
else
    invalid = value <= lowerLimit;
end
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || invalid
    if allowZero
        error('"%s" must be a non-negative numeric scalar.',name);
    else
        error('"%s" must be a positive numeric scalar.',name);
    end
end
end


function opts = parse_name_value_pairs(opts,varargin)
if mod(length(varargin),2) ~= 0
    error('Arguments must be provided as name-value pairs.');
end
for index = 1:2:length(varargin)
    name = lower(char(string(varargin{index})));
    if ~isfield(opts,name)
        error('Unknown option "%s".',varargin{index});
    end
    opts.(name) = varargin{index+1};
end
end
