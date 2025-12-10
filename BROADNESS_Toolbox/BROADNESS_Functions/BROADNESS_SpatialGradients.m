function [S_GRAD] = BROADNESS_SpatialGradients(BROADNESS, varargin)
%%
% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  SPATIAL GRADIENTS EMBEDDING ESTIMATION
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
%  This function computes the spatial gradients embedding of BROADNESS networks, 
%  clustering voxels based on participation in the different networks.
%
%  Specifically, it:
%   - Compute the thresholded activation patterns scaling the weigths
%     coefficient obtained with BROADNESS_NetworkEstimation function
%   - Clusters voxels in PC space using k-means (user-defined k range)
%   - Determines the optimal number of clusters using silhouette scores
%   - Saves cluster information, centroids, and NIFTI images (if path provided)
%   - Optionally generates 2D/3D scatterplots of voxels colored by cluster
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
%      - 'scatterplots'                    : Set to 'all' to plot cluster results for all k
%      - 'outpath'                         : If provided, saves NIFTI maps of clusters (default: [])
%      - 'mni_coords'                      : MNI coordinates (Nvoxels x 3) for 3D plotting in brain template
%                                            If empty, trying to read a default from files in 'External' function. This is in MNI space 8mm (LBPD order) 
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - S_GRAD                        : Structure with clustering results
%      - .idx                      : Table with cluster assignments (voxels × nclusters)
%      - .SUM                      : Table of within-cluster sums of distances
%      - .Centroids                : Cluster centroids for each k
%      - .optimalK                 : Optimal number of clusters (based on silhouette)
%
%
% ------------------------------------------------------------------------
%  NOTES:
% ------------------------------------------------------------------------
%
%  - The clustering is performed on z-scored, thresholded spatial maps.
%
%  - If 'outpath' is specified, the function saves NIFTI masks for each cluster
%    (only for the optimal k) using an 8mm MNI template.
%
%  - K-means clustering is repeated multiple times (Replicates = 100) for stability.
%
%  - The optimal number of clusters is determined as the mode of silhouette-based
%    evaluations repeated 10 times, to improve robustness.
%
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Chiara Malvaso, Mattia Rosso & Leonardo Bonetti 
%  chiara.malvaso@studio.unibo.it
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Department of Physics, University of Bologna
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 23/08/2025
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
    'scatterplots', [], ...     % [] (only optimal k), or 'all'
    'outpath', [], ...           % folder to save NIFTI masks (optional)
    'mni_coords', [] ...        % MNI coordinates for 3D plotting in brain template
);

% Parse name-value pairs
params = parse_name_value_pairs(params, varargin{:});

% Assign to readable internal names
selectedPCs              = params.principalcomps;
clusterRange             = params.nclusters;
plotMode                 = params.scatterplots;
activationThresh         = params.thresh;
savePath                 = params.outpath;
numSilhouetteRepeats     = params.evalclusters;
mni_coords               = params.mni_coords;

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

% Validate optional args
if ~isnumeric(selectedPCs) || ~isvector(selectedPCs)
    error('Please provide "principalcomps" as a numeric array or vector.');
end

if ~isvector(clusterRange)
    error('Please provide "nclusters" as a numeric vector.');
end

if isequal(selectedPCs, [1 2])
    disp('Computing the spatial gradient for 2 principal components (default).');
end

%% --------------- Compute thresholded activation patterns ----------------

disp('Computing Spatial Activation Patterns');

% Average across conditions (if 3D), to compute covariance on main effect
voxelTimeData = mean(voxelTimeData,4); %average across participants (if needed)
meanAcrossConditions = squeeze(mean(voxelTimeData(:,:,:), 3));  % voxels × time

% Covariance across voxels (rows are variables => transpose)
voxelCovariance = cov(meanAcrossConditions');

% Preallocate: thresholded activations used for clustering/plots
thresholdedActivations = zeros(nVoxels, length(selectedPCs));

% Loop over PCs up to the largest index in selectedPCs
for pcIdx = 1:selectedPCs(end)
    % Compute activation pattern for this PC (kept for completeness)
    % Note: not directly used downstream but retained to preserve original steps
    activationPatternForPC = activationWeights(:, pcIdx)' * voxelCovariance; %#ok<NASGU>

    % Determine threshold for voxel inclusion for this PC
    if isempty(activationThresh)
        pcThreshold = mean(abs(activationWeights(:, pcIdx))) + std(abs(activationWeights(:, pcIdx)));
    else
        warning('Using the input threshold for activation patterns. Ensure it suits your data.');
        pcThreshold = activationThresh;
    end

    % Apply threshold per voxel: keep original weight if above threshold, else 0
    for voxelIdx = 1:nVoxels
        if abs(activationWeights(voxelIdx, pcIdx)) > pcThreshold
            thresholdedActivations(voxelIdx, pcIdx) = activationWeights(voxelIdx, pcIdx);
        else
            thresholdedActivations(voxelIdx, pcIdx) = 0;
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

% Set the random seed for reproducibility
rng(42, 'twister');

for kIdx = 1:length(clusterRange)
    kVal = clusterRange(kIdx);

    % K-means with many replicates for stability
    % (Random initializations; same concept as original)
    [clusterLabels, centroids, sumD] = kmeans(zscoreActivations, kVal, 'Replicates', 100);

    % Store k and assignments
    clusterAssignmentsAll(1, kIdx)   = kVal;
    clusterAssignmentsAll(2:end,kIdx)= clusterLabels;

    % Store elbow metric
    withinClusterSums(kIdx,:) = [kVal, sum(sumD)];

    % Store centroids
    clusterCentroidsAll{kIdx} = centroids;
end

% Convert to user-friendly tables (columns labeled by k)
rawNames = strcat('Nclusters_', cellstr(num2str(clusterRange(:))));
varNamesByK = matlab.lang.makeValidName(rawNames);
S_GRAD.idx              = array2table(clusterAssignmentsAll(2:end,:), 'VariableNames', varNamesByK);
S_GRAD.SUM              = array2table(withinClusterSums(:,2)', 'VariableNames', varNamesByK);
S_GRAD.Centroids        = cell2table(clusterCentroidsAll, 'VariableNames', varNamesByK);

% -------- Elbow plot (sum of distances vs number of clusters) -----------
figure;

disp('Generating variance explained plot...');
plot(withinClusterSums(:,1), withinClusterSums(:,2), '-*', 'Color', [0.1882 0.4902 0.8118], 'LineWidth', 1.5);
hold on; grid minor; box on; set(gcf,'color','w');
xlabel('Number of clusters'); ylabel('SUM');
xlim([0, max(withinClusterSums(:,1)) + 1]);

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

S_GRAD.optimalK = optimalK;

%% ---------------- Prepare cluster-specific info for the optimal k -------

% Note: this section prepares cluster-wise voxel lists / tables only for
% the optimal solution; users can adapt to export as needed.

coordinates = load('MNI152_8mm_coord_dyi.mat'); % must contain coordinates.MNI8 (voxels × 3)
optimalCol  = find(clusterAssignmentsAll(1,:) == optimalK, 1, 'first');
clustersForOptimalK = clusterAssignmentsAll(2:end, optimalCol);   % voxel-wise labels 1..optimalK

% Build headers for a potential table per cluster: [VoxelIdx, X, Y, Z, PC1, PC2, ...]
tableHeaders = {'Voxel','X','Y','Z'};
for pcVal = selectedPCs
    tableHeaders{end+1} = ['PC' num2str(pcVal)];
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

S_GRAD.Clusters_info    = Clusters_info;

%% ----------------------- Scatter plots (PC space) -----------------------

% Default: plot only for optimal k. If 'all', iterate across all k in clusterRange.
if isempty(plotMode)
    if length(selectedPCs) == 2 || max(selectedPCs) == 3
        IDX = clustersForOptimalK;
        markerSize = 30;

        cmap = jet(optimalK) * 0.9;
        figure;
        for cl = 1:optimalK
            if length(selectedPCs) == 2
                scatter( ...
                    thresholdedActivations(IDX == cl, selectedPCs(1))', ...
                    thresholdedActivations(IDX == cl, selectedPCs(2))', ...
                    markerSize, 'MarkerFaceColor', cmap(cl,:), 'MarkerEdgeColor', cmap(cl,:) ...
                );
                hold on;
            else
                scatter3( ...
                    thresholdedActivations(IDX == cl, selectedPCs(1))', ...
                    thresholdedActivations(IDX == cl, selectedPCs(2))', ...
                    thresholdedActivations(IDX == cl, selectedPCs(3))', ...
                    markerSize, 'MarkerFaceColor', cmap(cl,:), 'MarkerEdgeColor', cmap(cl,:) ...
                );
                hold on;
            end
        end
        set(gcf,'color','w');
        grid minor; box on;
        colormap(cmap);
        cbar = colorbar;

        % Map colorbar limits to integer cluster labels (MATLAB-old compatible)
        if exist('clim', 'builtin')
            clim([0.5, optimalK + 0.5]);
        else
            caxis([0.5, optimalK + 0.5]);
        end
        cbar.Ticks = 1:optimalK;
        cbar.TickLabels = 1:optimalK;
        
        % ---- Axis labels ----
        ax = gca;
        xlabel(ax, ['Brain Network ' num2str(selectedPCs(1))]);
        ylabel(ax, ['Brain Network ' num2str(selectedPCs(2))]);
        if length(selectedPCs) >= 3
            zlabel(ax, ['Brain Network ' num2str(selectedPCs(3))]);
        end
        
    end
else
    if strcmp(plotMode, 'all')
        if ~(length(selectedPCs) == 2 || length(selectedPCs) == 3)
            warning('Scatter plots cannot be displayed for more than 3 dimensions.');
        else
            for kVal = clusterRange
                col = find(clusterAssignmentsAll(1,:) == kVal, 1, 'first');
                IDX = clusterAssignmentsAll(2:end, col);
                markerSize = 30;

                cmap = jet(kVal) * 0.9;
                figure;
                for cl = 1:kVal
                    if length(selectedPCs) == 2
                        scatter( ...
                            thresholdedActivations(IDX == cl, selectedPCs(1))', ...
                            thresholdedActivations(IDX == cl, selectedPCs(2))', ...
                            markerSize, 'MarkerFaceColor', cmap(cl,:), 'MarkerEdgeColor', cmap(cl,:) ...
                        );
                        hold on;
                    end

                    if length(selectedPCs) == 3
                        scatter3( ...
                            thresholdedActivations(IDX == cl, selectedPCs(1))', ...
                            thresholdedActivations(IDX == cl, selectedPCs(2))', ...
                            thresholdedActivations(IDX == cl, selectedPCs(3))', ...
                            markerSize, 'MarkerFaceColor', cmap(cl,:), 'MarkerEdgeColor', cmap(cl,:) ...
                        );
                        hold on;
                    end
                end
                set(gcf,'color','w');
                grid minor; box on;
                colormap(cmap);
                cbar = colorbar;
                
                % Map colorbar limits to integer cluster labels (MATLAB-old compatible)
                if exist('clim', 'builtin')
                    clim([0.5, kVal + 0.5]);   % sets range of color values
                else
                    caxis([0.5, kVal + 0.5]);
                end
                cbar.Ticks = 1:kVal;
                cbar.TickLabels = 1:kVal;
                % ---- New axis labels ----
                xlabel(['Brain Network ' num2str(selectedPCs(1))]);
                ylabel(['Brain Network ' num2str(selectedPCs(2))]);
                if length(selectedPCs) == 3
                    zlabel(['Brain Network ' num2str(selectedPCs(3))]);
                end
                
                title(['Cluster ' num2str(kVal)]);
            end
        end
    end
end

%% ------------------------- Save NIFTI images -----------------------

if ~isempty(savePath)
    disp('Saving NIFTI images only for the optimal number of clusters');
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
        save_nii(nii, [savePath '/BROADNESS_Output/BROADNESS_nifti/SpatialGradients_OptimalK_' num2str(optimalK) '_Cluster_' num2str(cl) '.nii']);
    end
end

%% ---------------- 3D brain per cluster (cluster membership only) --------

% Requires: 'mni_coords' as [Nvox x 3], same row order as data
if isempty(mni_coords)
    warning('Skipping 3D brain plots: provide ''mni_coords'' as [Nvox x 3].');
else
    MNIcoordsAll = mni_coords;
    if size(MNIcoordsAll,1) ~= nVoxels
        warning('Skipping 3D brain plots: mni_coords size mismatch (got %d rows, expected %d).', ...
                size(MNIcoordsAll,1), nVoxels);
    else
        disp('Generating 3D cluster maps (one brain per cluster, membership only)...');

        skipper       = 1;                % downsampling step
        scale_size    = 6;                % dot size (try 4–8 for visibility)
        templateFigFn = 'BrainTemplate_GT.fig';

        % actual cluster labels present
        labels   = clustersForOptimalK(:);
        uniqLabs = unique(labels(:)');    % e.g. could be [0 1 2] or [2 3 4]
        nLabs    = numel(uniqLabs);
        cmap     = jet(nLabs) * 0.9;      % one color per label

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
            if exist(templateFigFn,'file')
                fig = openfig(templateFigFn, 'new', 'visible');
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
                fig = figure('Color','w');
                ax  = axes('Parent', fig);
                hold(ax,'on');
            end

            % plot voxels for this cluster
            plot3(ax, coords(:,1), coords(:,2), coords(:,3), '.', ...
                  'Color', cmap(li,:), 'MarkerSize', scale_size);

            % styling / view
            axis(ax,'tight'); axis(ax,'equal'); axis(ax,'vis3d'); axis(ax,'off');
            rotate3d(ax,'on'); camlight(ax,'headlight'); lighting(ax,'gouraud');
            title(ax, sprintf('3D Cluster Map — OptimalK=%d — Cluster %g (n=%d voxels)', ...
                  optimalK, clLab, nnz(voxMask)), ...
                  'FontSize', 14, 'FontWeight', 'bold');
        end
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

