function [SPATIAL_CLUSTERING, FIGURES] = BROADNESS_SpatialActivationClustering(BROADNESS, varargin)
%%
% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  SPATIAL ACTIVATION PATTERN CLUSTERING
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
%  This function clusters the spatial activation patterns of BROADNESS
%  networks, grouping voxels according to their loading profiles across
%  the selected networks.
%
%  Specifically, it:
%   - Computes the thresholded activation patterns by scaling the weight
%     coefficient obtained with BROADNESS_NetworkEstimation function
%   - Clusters voxels in brain-network loading space using k-means
%     (user-defined k range)
%   - Determines the optimal number of clusters using silhouette scores
%   - Relabels clusters from the broadest network representation to the
%     least represented/background solution
%   - Saves cluster information, centroids, and NIFTI images (if path provided)
%   - Optionally generates 2D/3D scatterplots of voxels colored by cluster
%   - Generates one-dimensional voxel activation plots for the optimal clusters
%   - Optionally generates 3D plots in a brain template of the different clusters
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - BROADNESS                             : Structure outputted by BROADNESS_NetworkEstimation function.
%      - .OriginalData                     : 2D or 3D matrix (voxels × time × [conditions])
%      - .ActivationPatterns_BrainNetworks : 2D matrix (voxels × components)
%
%  - Optional arguments (name-value pairs):
%      - 'principalcomps'                  : Vector of PC indices to use (default: [1 2])
%      - 'nclusters'                       : Range of k-means clusters to test (default: 2:20)
%      - 'evalclusters'                    : Number of replications of the clustering analysis to identify ideal clustering solution using Silhouette method (default = 10) 
%      - 'thresh'                          : Threshold for including voxel activations (default: mean + std)
%      - 'representation_threshold'        : Minimum proportion of suprathreshold voxels required for a
%                                            network to be represented in a cluster (default: 0.50)
%      - 'scatterplots'                    : Set to 'all' to plot cluster results for all k
%      - 'outpath'                         : If provided, saves NIFTI maps and cluster activation plots (default: [])
%      - 'mni_coords'                      : MNI coordinates (Nvoxels x 3) for 3D plotting in brain template
%                                            If empty, trying to read a default from files in 'External' function. This is in MNI space 8mm (LBPD order)
%      - 'brainmarkersize'                 : Marker size for 3D cluster maps (default: 8)
%      - 'braincolorintensity'             : Multiplicative brightness of 3D cluster colors, between 0 and 1
%                                            (default: 0.75)
%      - 'figuremode'                      : 'off', 'show', 'save', or 'both' (default: 'off')
%      - 'figurelayout'                    : 'individual', 'summary', or 'both' (default: 'individual')
%      - 'figureoutpath'                   : Optional figure-only output folder. If empty, 'outpath'
%                                            is used when figures are saved (default: [])
%      - 'figureformats'                   : 'png', 'pdf', 'fig', or a cell array (default: {'png'})
%      - 'figureprefix'                    : Optional prefix for saved figure filenames
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - SPATIAL_CLUSTERING            : Structure with clustering results
%      - .idx                      : Table with cluster assignments (voxels × nclusters)
%      - .SUM                      : Table of within-cluster sums of distances
%      - .Centroids                : Cluster centroids for each k
%      - .optimalK                 : Optimal number of clusters (based on silhouette)
%      - .ClusterSummary           : Breadth, support, strength, direction, and voxel count for each
%                                    relabelled cluster in the optimal solution
%      - .ClusterSummaryAll        : Cluster summaries for every tested k
%      - .ClusterPoints_PC         : Per-cluster tables of voxel activations for the selected PCs
%      - .ClusterMinMax_PC         : Per-cluster minima and maxima for the selected PCs
%  - FIGURES                       : Visible figure handles and saved figure paths
%
%
% ------------------------------------------------------------------------
%  NOTES:
% ------------------------------------------------------------------------
%
%  - The clustering is performed on z-scored, thresholded spatial
%    activation patterns. Anatomical proximity between voxels is not used
%    by k-means.
%
%  - If 'outpath' is specified, the function saves NIFTI masks for each cluster
%    (only for the optimal k) using an 8mm MNI template.
%
%  - K-means clustering is repeated multiple times (Replicates = 100) for stability.
%
%  - K-means labels have no intrinsic order. Here they are relabelled after
%    clustering so Cluster 1 is the solution involving the broadest set of
%    selected networks and the background solution is placed last. This
%    changes only cluster names, not voxel membership or k-means results.
%
%  - Network representation is based on absolute suprathreshold activation.
%    Signed cluster means are retained separately in ClusterSummary.
%
%  - The bundled 1mm MNI152 template used for 3D cluster plots includes
%    the cerebellum and brainstem.
%
%  - The optimal number of clusters is determined as the mode of silhouette-based
%    evaluations repeated 10 times, to improve robustness.
%
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Chiara Malvaso, Mattia Rosso, Mathias Houe Andersen & Leonardo Bonetti 
%  chiara.malvaso@studio.unibo.it
%  mattia.rosso@clin.au.dk
%  mathias.houe.andersen@regionh.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Department of Physics, University of Bologna
%  Danish Research Centre for Magnetic Resonance, Copenhagen University Hospital
%  Faculty of Health and Medical Sciences, University of Copenhagen
%  Aarhus (DK), Copenhagen (DK), Oxford (UK), Bologna (Italy), Updated version 23/08/2025
%
% ========================================================================
%
%  This function uses the NIFTI Toolbox by Jimmy Shen:
%  Jimmy Shen (2025). Tools for NIFTI and ANALYZE image
%  (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
%  MATLAB Central File Exchange. Retrieved July 05, 2025.

%%





%% ----------------------------- Parse inputs -----------------------------

disp('Checking inputs');

% Defaults
params = struct( ...
    'principalcomps', 1:2, ...  % PCs to use in the embedding/plots
    'nclusters', 2:20, ...      % range of k for k-means
    'evalclusters',10, ...      % repetitions for clustering analysis to identify ideal clustering solution
    'thresh', [], ...           % per-PC abs(weight) threshold; default mean+std
    'representation_threshold', 0.50, ... % minimum within-cluster support for a represented network
    'scatterplots', [], ...     % [] (only optimal k), or 'all'
    'outpath', [], ...           % folder to save NIFTI masks (optional)
    'mni_coords', [], ...       % MNI coordinates for 3D plotting in brain template
    'brainmarkersize', 8, ...   % marker size for 3D cluster maps
    'braincolorintensity', 0.75, ... % brightness multiplier for 3D cluster colors
    'figuremode', 'off', ...
    'figurelayout', 'individual', ...
    'figureoutpath', [], ...
    'figureformats', {{'png'}}, ...
    'figureprefix', '' ...
);

% Parse name-value pairs
params = parse_name_value_pairs(params, varargin{:});

% Assign to readable internal names
selectedPCs              = params.principalcomps;
clusterRange             = params.nclusters;
plotMode                 = params.scatterplots;
activationThresh         = params.thresh;
representationThreshold  = params.representation_threshold;
savePath                 = params.outpath;
numSilhouetteRepeats     = params.evalclusters;
mni_coords               = params.mni_coords;
brainMarkerSize          = params.brainmarkersize;
brainColorIntensity      = params.braincolorintensity;
figureOutputPath         = params.figureoutpath;
if isempty(figureOutputPath), figureOutputPath = savePath; end
figureSettings = BROADNESS_FigureSettings(params.figuremode, params.figurelayout, ...
    figureOutputPath, params.figureformats, params.figureprefix, 'SpatialActivationClustering');
figureHandles = gobjects(0);
figureFiles = {};

% Validate required fields
if isfield(BROADNESS, 'OriginalData')
    voxelTimeData = BROADNESS.OriginalData;  % voxels × time × [conditions]
else
    error('Invalid input: BROADNESS.OriginalData is required.');
end

% if ndims(voxelTimeData) ~= 3 && ~ismatrix(voxelTimeData)
%     error('"OriginalData" must be a 2D or 3D matrix: (voxels × time × [conditions]).');
% end

nVoxels = size(voxelTimeData, 1);

if isfield(BROADNESS, 'ActivationPatterns_BrainNetworks')
    activationWeights = BROADNESS.ActivationPatterns_BrainNetworks; % voxels × components
else
    error('Invalid input: BROADNESS.ActivationPatterns_BrainNetworks is required.');
end

if ~ismatrix(activationWeights)
    error('"ActivationPatterns_BrainNetworks" must be a 2D matrix (voxels × components).');
end
if size(activationWeights,1) ~= nVoxels
    error('The number of voxels in "OriginalData" and "ActivationPatterns_BrainNetworks" must match.');
end

% Validate optional args
if ~isnumeric(selectedPCs) || ~isvector(selectedPCs) || isempty(selectedPCs)
    error('Please provide "principalcomps" as a non-empty numeric vector.');
end
selectedPCs = selectedPCs(:)'; % use a consistent row-vector representation
if any(~isfinite(selectedPCs)) || any(selectedPCs < 1) || any(fix(selectedPCs) ~= selectedPCs)
    error('The values in "principalcomps" must be positive integer component indices.');
end
if any(selectedPCs > size(activationWeights,2))
    error('The requested "principalcomps" exceed the available brain-network components.');
end
if length(unique(selectedPCs)) ~= length(selectedPCs)
    error('The values in "principalcomps" must be unique.');
end

if ~isvector(clusterRange)
    error('Please provide "nclusters" as a numeric vector.');
end
if ~isnumeric(representationThreshold) || ~isscalar(representationThreshold) || ...
        ~isfinite(representationThreshold) || representationThreshold <= 0 || representationThreshold > 1
    error('"representation_threshold" must be a numeric scalar greater than 0 and no larger than 1.');
end
if ~isnumeric(brainMarkerSize) || ~isscalar(brainMarkerSize) || ...
        ~isfinite(brainMarkerSize) || brainMarkerSize <= 0
    error('"brainmarkersize" must be a positive numeric scalar.');
end
if ~isnumeric(brainColorIntensity) || ~isscalar(brainColorIntensity) || ...
        ~isfinite(brainColorIntensity) || brainColorIntensity <= 0 || brainColorIntensity > 1
    error('"braincolorintensity" must be a numeric scalar greater than 0 and no larger than 1.');
end

if isequal(selectedPCs, [1 2])
    disp('Computing spatial activation clustering for 2 principal components (default).');
end

%% --------------- Compute thresholded activation patterns ----------------

disp('Computing Spatial Activation Patterns');

% Preallocate: thresholded activations used for clustering/plots
thresholdedActivations = zeros(nVoxels, length(selectedPCs));
activationThresholds = zeros(1, length(selectedPCs));

% Loop over the selected PCs, storing them in the requested order
for pcCol = 1:length(selectedPCs)
    pcIdx = selectedPCs(pcCol);
    % Determine threshold for voxel inclusion for this PC
    if isempty(activationThresh)
        pcThreshold = mean(abs(activationWeights(:, pcIdx))) + std(abs(activationWeights(:, pcIdx)));
    else
        warning('Using the input threshold for activation patterns. Ensure it suits your data.');
        pcThreshold = activationThresh;
    end
    activationThresholds(pcCol) = pcThreshold;

    % Apply threshold per voxel: keep original weight if above threshold, else 0
    for voxelIdx = 1:nVoxels
        if abs(activationWeights(voxelIdx, pcIdx)) > pcThreshold
            thresholdedActivations(voxelIdx, pcCol) = activationWeights(voxelIdx, pcIdx);
        else
            thresholdedActivations(voxelIdx, pcCol) = 0;
        end
    end
end

%% --------------------------- K-means clustering -------------------------

% Z-score across voxels so PCs are comparable in scale
zscoreActivations = (thresholdedActivations - mean(thresholdedActivations)) ./ std(thresholdedActivations);

% Prepare storage across the tested k values
clusterAssignmentsAll = zeros(nVoxels + 1, length(clusterRange)); % first row stores k itself
withinClusterSums     = zeros(length(clusterRange), 2);            % [k, sum(sumD)]
clusterCentroidsAll   = cell(1, length(clusterRange));             % centroids per k
clusterSummariesAll   = cell(1, length(clusterRange));             % breadth-based summary per k

% Set the random seed for reproducibility
rng(42, 'twister');

for kIdx = 1:length(clusterRange)
    kVal = clusterRange(kIdx);

    % K-means with many replicates for stability
    % (Random initializations; same concept as original)
    [clusterLabels, centroids, sumD] = kmeans(zscoreActivations, kVal, 'Replicates', 100);

    % K-means cluster numbers are arbitrary. Relabel each solution so the
    % broadest network representations appear first and background last.
    [clusterLabels, centroids, sumD, clusterSummary] = order_clusters_by_network_breadth( ...
        clusterLabels, centroids, sumD, thresholdedActivations, ...
        selectedPCs, representationThreshold);

    % Store k and assignments
    clusterAssignmentsAll(1, kIdx)   = kVal;
    clusterAssignmentsAll(2:end,kIdx)= clusterLabels;

    % Store elbow metric
    withinClusterSums(kIdx,:) = [kVal, sum(sumD)];

    % Store centroids
    clusterCentroidsAll{kIdx} = centroids;
    clusterSummariesAll{kIdx} = clusterSummary;
end

% Convert to user-friendly tables (columns labeled by k)
rawNames = strcat('Nclusters_', cellstr(num2str(clusterRange(:))));
varNamesByK = matlab.lang.makeValidName(rawNames);
SPATIAL_CLUSTERING.idx       = array2table(clusterAssignmentsAll(2:end,:), 'VariableNames', varNamesByK);
SPATIAL_CLUSTERING.SUM       = array2table(withinClusterSums(:,2)', 'VariableNames', varNamesByK);
SPATIAL_CLUSTERING.Centroids = cell2table(clusterCentroidsAll, 'VariableNames', varNamesByK);
SPATIAL_CLUSTERING.ClusterSummaryAll = cell2table(clusterSummariesAll, 'VariableNames', varNamesByK);
SPATIAL_CLUSTERING.RepresentationThreshold = representationThreshold;
SPATIAL_CLUSTERING.ActivationThresholds = array2table(activationThresholds, ...
    'VariableNames', matlab.lang.makeValidName(strcat('BN', cellstr(num2str(selectedPCs(:))))));
SPATIAL_CLUSTERING.SelectedComponents = selectedPCs;

% -------- Elbow plot (sum of distances vs number of clusters) -----------
if ~strcmp(figureSettings.Mode, 'off') && figureSettings.MakeIndividual
    disp('Generating clustering elbow plot...');
    fig = figure('Visible', figureSettings.Visible, 'Color', 'w');
    plot_clustering_elbow(gca, withinClusterSums);
    figureFiles = [figureFiles; BROADNESS_FinalizeFigure(fig, figureSettings, ...
        'ClusteringElbow')];
    if figureSettings.Show, figureHandles(end+1) = fig; end
end

%% -------------------- Determine optimal number of clusters --------------

% Repeat silhouette evaluation to stabilize selection; choose the mode
disp('Computing the optimal number of clusters with evalclusters');

if length(clusterRange) > 1
    optimalKAcrossRepeats = zeros(numSilhouetteRepeats, 1);
    for rep = 1:numSilhouetteRepeats
        rng(rep, 'twister'); % Different seed per repeat, but reproducible overall
        evaObj = evalclusters(zscoreActivations, 'kmeans', 'silhouette', 'KList', clusterRange);
        optimalKAcrossRepeats(rep) = evaObj.OptimalK;
        disp(['Computing the optimal number of clusters with evalclusters - Repetition ' num2str(rep) ' / ' num2str(numSilhouetteRepeats)])
    end
    optimalK = mode(optimalKAcrossRepeats);
else
    disp('Only one "nclusters" value provided — taking it as optimal.');
    optimalK = clusterRange;
end

SPATIAL_CLUSTERING.optimalK = optimalK;

%% ---------------- Prepare cluster-specific info for the optimal k -------

% Note: this section prepares cluster-wise voxel lists / tables only for
% the optimal solution; users can adapt to export as needed.

coordinates = load('MNI152_8mm_coord_dyi.mat'); % must contain coordinates.MNI8 (voxels × 3)
optimalCol  = find(clusterAssignmentsAll(1,:) == optimalK, 1, 'first');
clustersForOptimalK = clusterAssignmentsAll(2:end, optimalCol);   % voxel-wise labels 1..optimalK
SPATIAL_CLUSTERING.ClusterSummary = clusterSummariesAll{optimalCol};

% Build headers for a potential table per cluster: [VoxelIdx, X, Y, Z, PC1, PC2, ...]
tableHeaders = cell(1, 4 + length(selectedPCs));
tableHeaders(1:4) = {'Voxel','X','Y','Z'};
for pcCol = 1:length(selectedPCs)
    tableHeaders{4 + pcCol} = ['PC' num2str(selectedPCs(pcCol))];
end

Clusters_info = cell(optimalK,1);
% For each cluster, assemble a table of its voxels
for cl = 1:optimalK
    % Binary membership vector for this cluster
    clusterMaskBinary = (clustersForOptimalK == cl);

    % Collect activations and coordinates for voxels in this cluster
    voxelIdxList = find(clusterMaskBinary);
    clusterActivations = thresholdedActivations(voxelIdxList, :);
    if isempty(mni_coords)
        MNIcoords = coordinates.MNI8(voxelIdxList, :);
    else
        MNIcoords = mni_coords(voxelIdxList, :);
    end

    % Compose table data: a running index within the cluster, coords, activations
    excel_data = [(1:size(clusterActivations,1))', MNIcoords, clusterActivations];
    tbl = array2table(excel_data, 'VariableNames', tableHeaders);
    Clusters_info{cl} = tbl; % Store information
end

SPATIAL_CLUSTERING.Clusters_info = Clusters_info;

%% ------------ One-dimensional activation plots by cluster -------------

% Store the voxel values used in each cluster plot. The columns of
% thresholdedActivations already follow the order requested in selectedPCs.
nSelectedPCs = length(selectedPCs);
pcVariableNames = cell(1, nSelectedPCs);
minMaxVariableNames = cell(1, 1 + 2 * nSelectedPCs);
minMaxVariableNames{1} = 'Cluster';
for pcCol = 1:nSelectedPCs
    pcVariableNames{pcCol} = ['PC' num2str(selectedPCs(pcCol))];
    minMaxVariableNames{2 * pcCol} = ['PC' num2str(selectedPCs(pcCol)) '_min'];
    minMaxVariableNames{2 * pcCol + 1} = ['PC' num2str(selectedPCs(pcCol)) '_max'];
end

ClusterPoints_PC = cell(optimalK, 1);
clusterMinMaxValues = zeros(optimalK, 1 + 2 * nSelectedPCs);
clusterMinMaxValues(:,1) = (1:optimalK)';

for cl = 1:optimalK
    clusterValues = thresholdedActivations(clustersForOptimalK == cl, :);
    ClusterPoints_PC{cl} = array2table(clusterValues, ...
        'VariableNames', pcVariableNames);

    for pcCol = 1:nSelectedPCs
        clusterMinMaxValues(cl, 2 * pcCol) = min(clusterValues(:,pcCol));
        clusterMinMaxValues(cl, 2 * pcCol + 1) = max(clusterValues(:,pcCol));
    end
end

ClusterMinMax_PC = array2table(clusterMinMaxValues, ...
    'VariableNames', minMaxVariableNames);
SPATIAL_CLUSTERING.ClusterPoints_PC = ClusterPoints_PC;
SPATIAL_CLUSTERING.ClusterMinMax_PC = ClusterMinMax_PC;

% Use a common symmetric x-axis so the optimal clusters are comparable.
finiteActivationValues = thresholdedActivations(isfinite(thresholdedActivations));
if isempty(finiteActivationValues)
    maxAbsoluteActivation = 1;
else
    maxAbsoluteActivation = max(abs(finiteActivationValues));
    if maxAbsoluteActivation == 0
        maxAbsoluteActivation = 1;
    end
end
plotLimits = [-1 1] * maxAbsoluteActivation * 1.05;

clusterColors = jet(optimalK) * 0.9;
laneSeparation = 1;
lanePositions = (nSelectedPCs-1:-1:0) * laneSeparation;
pointJitter = 0.28 * laneSeparation;
baselineColor = [0.65 0.65 0.65];

if ~strcmp(figureSettings.Mode, 'off')
disp('Generating one-dimensional activation plots for the optimal clusters');
activationSummaryFigure = [];
activationSummaryLayout = [];
if figureSettings.MakeSummary
    activationSummaryFigure = figure('Visible', figureSettings.Visible, ...
        'Color', 'w', 'Position', [100 100 1200 750]);
    activationSummaryLayout = tiledlayout(activationSummaryFigure, 'flow', ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    title(activationSummaryLayout, 'Optimal Clusters — Network Activation Profiles');
end
for cl = 1:optimalK
    clusterTable = ClusterPoints_PC{cl};
    clusterVisibility = figureSettings.Visible;
    if ~figureSettings.MakeIndividual, clusterVisibility = 'off'; end
    clusterFigure = figure('Color', 'w', ...
        'Visible', clusterVisibility, ...
        'Name', ['Cluster ' num2str(cl) ' activation values'], ...
        'NumberTitle', 'off');
    clusterAxes = axes('Parent', clusterFigure);
    hold(clusterAxes, 'on');

    for pcCol = 1:nSelectedPCs
        laneY = lanePositions(pcCol);
        plot(clusterAxes, plotLimits, [laneY laneY], '-', ...
            'Color', baselineColor, 'LineWidth', 1);
        plot(clusterAxes, [plotLimits(1) plotLimits(1)], ...
            laneY + [-0.15 0.15] * laneSeparation, '-', ...
            'Color', baselineColor, 'LineWidth', 1);
        plot(clusterAxes, [plotLimits(2) plotLimits(2)], ...
            laneY + [-0.15 0.15] * laneSeparation, '-', ...
            'Color', baselineColor, 'LineWidth', 1);

        text(clusterAxes, plotLimits(1) - 0.02 * diff(plotLimits), laneY, ...
            ['BN' num2str(selectedPCs(pcCol))], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

        activationValues = clusterTable{:,pcCol};
        activationValues = activationValues(isfinite(activationValues));
        % Exact zeros represent values removed by thresholding. Retain
        % them in the output tables, but omit them from the plot so they
        % do not form an artificial vertical column at zero.
        activationValues = activationValues(activationValues ~= 0);
        if ~isempty(activationValues)
            % Continuous deterministic jitter separates overlapping points
            % without creating artificial rows or changing MATLAB's
            % random-number state.
            jitterOrder = mod((1:length(activationValues))' * ...
                0.618033988749895, 1);
            jitterValues = (2 * jitterOrder - 1) * pointJitter;
            scatter(clusterAxes, activationValues, laneY + jitterValues, 18, ...
                'MarkerFaceColor', clusterColors(cl,:), ...
                'MarkerEdgeColor', 'none');
        end
    end

    xlim(clusterAxes, plotLimits);
    ylim(clusterAxes, [min(lanePositions)-0.4 max(lanePositions)+0.4]);
    xticks(clusterAxes, [plotLimits(1) 0 plotLimits(2)]);
    set(clusterAxes, 'YColor', 'none');
    box(clusterAxes, 'off');
    grid(clusterAxes, 'off');
    title(clusterAxes, ['Cluster ' num2str(cl) ' — ' ...
        SPATIAL_CLUSTERING.ClusterSummary.NetworkCombination{cl}]);

    if figureSettings.MakeSummary
        summaryAxes = nexttile(activationSummaryLayout);
        copy_axes_content(clusterAxes, summaryAxes, ['Cluster ' num2str(cl)]);
        set(summaryAxes, 'YColor', 'none'); box(summaryAxes, 'off'); grid(summaryAxes, 'off');
    end
    if figureSettings.MakeIndividual
        figureFiles = [figureFiles; BROADNESS_FinalizeFigure(clusterFigure, ...
            figureSettings, ['OptimalK_' num2str(optimalK,'%02d') ...
            '_Cluster_' num2str(cl,'%02d') '_ActivationProfile'])]; %#ok<AGROW>
        if figureSettings.Show, figureHandles(end+1) = clusterFigure; end %#ok<AGROW>
    else
        close(clusterFigure)
    end
end
if figureSettings.MakeSummary
    figureFiles = [figureFiles; BROADNESS_FinalizeFigure(activationSummaryFigure, ...
        figureSettings, ['OptimalK_' num2str(optimalK,'%02d') ...
        '_ActivationProfiles_Summary'])];
    if figureSettings.Show, figureHandles(end+1) = activationSummaryFigure; end
end
end

%% ----------------------- Scatter plots (PC space) -----------------------

% Default: plot only the optimal k. 'scatterplots','all' produces an
% individual embedding for every tested k.
if ~strcmp(figureSettings.Mode, 'off')
    if ~(length(selectedPCs) == 2 || length(selectedPCs) == 3)
        warning('Cluster embeddings can be displayed only for 2 or 3 dimensions.');
    else
        if figureSettings.MakeIndividual
            if isempty(plotMode)
                kValuesToPlot = optimalK;
            elseif strcmpi(plotMode, 'all')
                kValuesToPlot = clusterRange;
            else
                kValuesToPlot = optimalK;
            end
            for kVal = kValuesToPlot
                col = find(clusterAssignmentsAll(1,:) == kVal, 1, 'first');
                fig = figure('Visible', figureSettings.Visible, 'Color', 'w');
                plot_cluster_embedding(gca, thresholdedActivations, ...
                    clusterAssignmentsAll(2:end,col), selectedPCs, kVal);
                figureFiles = [figureFiles; BROADNESS_FinalizeFigure(fig, ...
                    figureSettings, ['ClusterEmbedding_K_' num2str(kVal,'%02d')])]; %#ok<AGROW>
                if figureSettings.Show, figureHandles(end+1) = fig; end %#ok<AGROW>
            end
        end

        if figureSettings.MakeSummary
            fig = figure('Visible', figureSettings.Visible, 'Color', 'w', ...
                'Position', [100 100 1100 480]);
            layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            title(layout, ['Spatial Activation Clustering — Optimal K = ' num2str(optimalK)]);
            plot_clustering_elbow(nexttile(layout), withinClusterSums);
            plot_cluster_embedding(nexttile(layout), thresholdedActivations, ...
                clustersForOptimalK, selectedPCs, optimalK);
            figureFiles = [figureFiles; BROADNESS_FinalizeFigure(fig, ...
                figureSettings, ['OptimalK_' num2str(optimalK,'%02d') ...
                '_ClusteringSummary'])];
            if figureSettings.Show, figureHandles(end+1) = fig; end
        end
    end
end

%% ------------------------- Save NIFTI images -----------------------

if ~isempty(savePath)
    disp('Saving NIFTI images only for the optimal number of clusters');
    niftiPath = fullfile(savePath, 'BROADNESS_Output', 'BROADNESS_nifti');
    if ~exist(niftiPath, 'dir')
        mkdir(niftiPath);
    end
    maskNii = load_nii('MNI152_8mm_brain_diy.nii.gz'); % 8mm brain mask indexed by voxel id
    optimalCol  = find(clusterAssignmentsAll(1,:) == optimalK, 1, 'first');
    clustersForOptimalK = clusterAssignmentsAll(2:end, optimalCol); % voxel labels 1..optimalK

    % Build a binary mask per cluster and write as NIFTI
    volSize = size(maskNii.img);
    for cl = 1:optimalK
        % Binary membership vector for this cluster
        clusterMaskBinary = (clustersForOptimalK == cl);

        % Initialize empty 3D volume (no time dimension)
        outVol = zeros(volSize(1), volSize(2), volSize(3), 1);

        % For each voxel id, put its binary label into the voxel position
        for voxelId = 1:length(clusterMaskBinary)
            idxInMask = find(maskNii.img == voxelId);     % linear indices for this voxelId
            [i1,i2,i3] = ind2sub(volSize(1:3), idxInMask);
            outVol(i1,i2,i3,:) = clusterMaskBinary(voxelId,:);
        end

        % Make and save NIFTI (8 mm resolution preserved)
        nii = make_nii(outVol, [8 8 8]);
        nii.img = outVol;
        nii.hdr.hist = maskNii.hdr.hist; % copy header info
        disp(['Saving NIFTI image - cluster ' num2str(cl)]);
        save_nii(nii, fullfile(niftiPath, ...
            ['SpatialActivationClustering_OptimalK_' num2str(optimalK) ...
            '_Cluster_' num2str(cl) '.nii']));
    end
end

%% ---------------- 3D brain per cluster (cluster membership only) --------

% Requires: 'mni_coords' as [Nvox x 3], same row order as data
if strcmp(figureSettings.Mode, 'off')
    % Numerical clustering and NIFTI outputs remain available without plots.
elseif isempty(mni_coords)
    warning('Skipping 3D brain plots: provide ''mni_coords'' as [Nvox x 3].');
else
    MNIcoordsAll = mni_coords;
    if size(MNIcoordsAll,1) ~= nVoxels
        warning('Skipping 3D brain plots: mni_coords size mismatch (got %d rows, expected %d).', ...
                size(MNIcoordsAll,1), nVoxels);
    else
        disp('Generating 3D cluster maps (one brain per cluster, membership only)...');

        skipper       = 1;                % downsampling step
        scale_size    = brainMarkerSize;  % user-adjustable dot size
        templateFigFn = 'BrainTemplate_MNI152_1mm_FullBrain.fig';

        % actual cluster labels present
        labels   = clustersForOptimalK(:);
        uniqLabs = unique(labels(:)');    % e.g. could be [0 1 2] or [2 3 4]
        nLabs    = numel(uniqLabs);
        cmap     = jet(nLabs) * brainColorIntensity; % one color per label

        brainSummaryFigure = [];
        brainSummaryLayout = [];
        if figureSettings.MakeSummary
            brainSummaryFigure = figure('Visible', figureSettings.Visible, ...
                'Color', 'w', 'Position', [100 100 1200 750]);
            brainSummaryLayout = tiledlayout(brainSummaryFigure, 'flow', ...
                'TileSpacing', 'compact', 'Padding', 'compact');
            title(brainSummaryLayout, '3D Cluster Maps');
        end

        for li = 1:nLabs
            clLab   = uniqLabs(li);
            voxMask = (labels == clLab);

            if ~any(voxMask)
                warning('Cluster label %g: no voxels to plot.', clLab);
                continue;
            end

            coords = MNIcoordsAll(voxMask, :);
            coords = coords(1:skipper:end, :);

            % open template brain or new figure
            brainVisibility = figureSettings.Visible;
            if ~figureSettings.MakeIndividual, brainVisibility = 'off'; end
            if exist(templateFigFn,'file')
                if strcmp(brainVisibility, 'on')
                    openFigureVisibility = 'visible';
                else
                    openFigureVisibility = 'invisible';
                end
                fig = openfig(templateFigFn, 'new', openFigureVisibility);
                ax  = findobj(fig, 'Type','axes');
                if isempty(ax), ax = axes('Parent', fig); end
                ax = ax(1);
                set(fig,'Renderer','opengl');
                hold(ax,'on');

                % make template brain transparent if it's a patch object
%                 brainPatch = findobj(ax, 'Type','patch');
%                 if ~isempty(brainPatch)
%                     set(brainPatch, 'FaceAlpha', 0.15, 'EdgeColor','none');
%                 end
            else
                fig = figure('Color','w', 'Visible', brainVisibility);
                ax  = axes('Parent', fig);
                hold(ax,'on');
            end

            % plot voxels for this cluster
            plot3(ax, coords(:,1), coords(:,2), coords(:,3), '.', ...
                  'Color', cmap(li,:), 'MarkerSize', scale_size);

            % styling / view
            axis(ax,'tight'); axis(ax,'equal'); axis(ax,'vis3d'); axis(ax,'off');
            rotate3d(ax,'on'); camlight(ax,'headlight'); lighting(ax,'gouraud');
            title(ax, sprintf('3D Cluster Map — OptimalK=%d — Cluster %g: %s (n=%d voxels)', ...
                  optimalK, clLab, SPATIAL_CLUSTERING.ClusterSummary.NetworkCombination{clLab}, nnz(voxMask)), ...
                  'FontSize', 14, 'FontWeight', 'bold');

            if figureSettings.MakeSummary
                summaryAxes = nexttile(brainSummaryLayout, li);
                copy_axes_content(ax, summaryAxes, ['Cluster ' num2str(clLab)]);
                axis(summaryAxes, 'equal'); axis(summaryAxes, 'vis3d'); axis(summaryAxes, 'off');
                view(summaryAxes, 3); camlight(summaryAxes, 'headlight'); lighting(summaryAxes, 'gouraud');
                clusterTitle = title(summaryAxes, ['Cluster ' num2str(clLab)], ...
                    'FontWeight', 'bold');
                set(clusterTitle, 'Units', 'normalized', 'Position', [0.5 0.94 0]);
            end
            if figureSettings.MakeIndividual
                figureFiles = [figureFiles; BROADNESS_FinalizeFigure(fig, ...
                    figureSettings, ['OptimalK_' num2str(optimalK,'%02d') ...
                    '_Cluster_' num2str(clLab,'%02d') '_BrainMap'])]; %#ok<AGROW>
                if figureSettings.Show, figureHandles(end+1) = fig; end %#ok<AGROW>
            else
                close(fig)
            end
        end
        if figureSettings.MakeSummary
            figureFiles = [figureFiles; BROADNESS_FinalizeFigure(brainSummaryFigure, ...
                figureSettings, ['OptimalK_' num2str(optimalK,'%02d') ...
                '_BrainMaps_Summary'])];
            if figureSettings.Show, figureHandles(end+1) = brainSummaryFigure; end
        end
    end
end

FIGURES.Handles = figureHandles;
FIGURES.Files = figureFiles;
SPATIAL_CLUSTERING.Figures.Files = figureFiles;

%% ------------------------ Helper: parse name/values ---------------------

function plot_clustering_elbow(ax, withinClusterSums)
plot(ax, withinClusterSums(:,1), withinClusterSums(:,2), '-*', ...
    'Color', [0.1882 0.4902 0.8118], 'LineWidth', 1.5, 'MarkerSize', 7);
grid(ax, 'minor'); box(ax, 'on');
xlabel(ax, 'Number of clusters'); ylabel(ax, 'Within-cluster sum of distances');
xlim(ax, [0, max(withinClusterSums(:,1)) + 1]); title(ax, 'Clustering elbow');
end

function plot_cluster_embedding(ax, activations, labels, selectedPCs, numberClusters)
hold(ax, 'on');
colors = jet(numberClusters) * 0.9;
for clusteri = 1:numberClusters
    if length(selectedPCs) == 2
        scatter(ax, activations(labels == clusteri,1), activations(labels == clusteri,2), ...
            24, 'MarkerFaceColor', colors(clusteri,:), 'MarkerEdgeColor', 'none');
    else
        scatter3(ax, activations(labels == clusteri,1), activations(labels == clusteri,2), ...
            activations(labels == clusteri,3), 24, 'MarkerFaceColor', colors(clusteri,:), ...
            'MarkerEdgeColor', 'none');
    end
end
grid(ax, 'minor'); box(ax, 'on'); colormap(ax, colors);
cbar = colorbar(ax); caxis(ax, [0.5 numberClusters+0.5]);
cbar.Ticks = 1:numberClusters; cbar.TickLabels = 1:numberClusters;
xlabel(ax, ['Brain Network ' num2str(selectedPCs(1))]);
ylabel(ax, ['Brain Network ' num2str(selectedPCs(2))]);
if length(selectedPCs) == 3
    zlabel(ax, ['Brain Network ' num2str(selectedPCs(3))]); view(ax, 3);
end
title(ax, ['Cluster embedding — K = ' num2str(numberClusters)]);
end

function copy_axes_content(sourceAxes, targetAxes, targetTitle)
copyobj(sourceAxes.Children, targetAxes);
set(targetAxes, 'XLim', sourceAxes.XLim, 'YLim', sourceAxes.YLim);
if isprop(sourceAxes, 'ZLim'), set(targetAxes, 'ZLim', sourceAxes.ZLim); end
xlabel(targetAxes, sourceAxes.XLabel.String);
ylabel(targetAxes, sourceAxes.YLabel.String);
if ~isempty(sourceAxes.ZLabel.String), zlabel(targetAxes, sourceAxes.ZLabel.String); end
title(targetAxes, targetTitle); box(targetAxes, sourceAxes.Box); grid(targetAxes, 'off');
end

function [newLabels, newCentroids, newSumD, summary] = order_clusters_by_network_breadth( ...
        oldLabels, oldCentroids, oldSumD, thresholdedActivations, ...
        selectedPCs, representationThreshold)

    nClusters = size(oldCentroids,1);
    nNetworks = length(selectedPCs);
    support = zeros(nClusters,nNetworks);
    normalizedStrength = zeros(nClusters,nNetworks);
    signedMean = zeros(nClusters,nNetworks);
    voxelCount = zeros(nClusters,1);
    firstVoxel = zeros(nClusters,1);

    networkScale = max(abs(thresholdedActivations),[],1);
    networkScale(networkScale == 0) = 1;

    for cluster = 1:nClusters
        clusterMask = oldLabels == cluster;
        clusterValues = thresholdedActivations(clusterMask,:);
        voxelCount(cluster) = sum(clusterMask);
        firstVoxel(cluster) = find(clusterMask,1,'first');
        support(cluster,:) = mean(clusterValues ~= 0,1);
        normalizedStrength(cluster,:) = mean(abs(clusterValues),1) ./ networkScale;
        signedMean(cluster,:) = mean(clusterValues,1);
    end

    representedNetworks = support >= representationThreshold;
    networkBreadth = sum(representedNetworks,2);
    totalSupport = sum(support,2);
    totalStrength = sum(normalizedStrength,2);

    % Sort descending by breadth, total support, strength, and individual
    % network representation. The first voxel is a deterministic final
    % tie-breaker that is independent of the arbitrary original label.
    sortValues = [networkBreadth totalSupport totalStrength ...
        representedNetworks normalizedStrength -firstVoxel];
    [~,newToOld] = sortrows(sortValues, -(1:size(sortValues,2)));

    oldToNew = zeros(nClusters,1);
    oldToNew(newToOld) = 1:nClusters;
    newLabels = oldToNew(oldLabels);
    newCentroids = oldCentroids(newToOld,:);
    newSumD = oldSumD(newToOld,:);

    networkCombination = cell(nClusters,1);
    for newCluster = 1:nClusters
        oldCluster = newToOld(newCluster);
        representedPCs = selectedPCs(representedNetworks(oldCluster,:));
        if isempty(representedPCs)
            networkCombination{newCluster} = 'None/background';
        else
            networkNames = arrayfun(@(pc) ['BN' num2str(pc)], ...
                representedPCs, 'UniformOutput', false);
            networkCombination{newCluster} = strjoin(networkNames,' + ');
        end
    end

    summary = table((1:nClusters)', newToOld(:), voxelCount(newToOld), ...
        networkBreadth(newToOld), totalSupport(newToOld), ...
        totalStrength(newToOld), networkCombination, ...
        'VariableNames', {'Cluster','OriginalCluster','VoxelCount', ...
        'RepresentedNetworks','TotalSupport','TotalNormalizedStrength', ...
        'NetworkCombination'});

    for network = 1:nNetworks
        networkName = ['BN' num2str(selectedPCs(network))];
        summary.([networkName '_Represented']) = representedNetworks(newToOld,network);
        summary.([networkName '_Support']) = support(newToOld,network);
        summary.([networkName '_NormalizedStrength']) = normalizedStrength(newToOld,network);
        summary.([networkName '_SignedMean']) = signedMean(newToOld,network);
    end
end

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

