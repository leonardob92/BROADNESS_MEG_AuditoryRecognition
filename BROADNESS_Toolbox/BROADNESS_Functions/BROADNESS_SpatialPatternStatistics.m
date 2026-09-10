function [SPATIAL_STATS, FIGURES] = BROADNESS_SpatialPatternStatistics(BROADNESS, varargin)
%%
% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  SPATIAL ACTIVATION PATTERN STATISTICS
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
%  This function compares condition-specific spatial activation patterns
%  across participants. Statistical inference is performed on unthresholded
%  maps; activation thresholds used for visualization must not be applied
%  before calling this function.
%
%  Specifically, it:
%   - Performs paired condition comparisons at every brain source
%   - Corrects the spatial tests with cluster-mass permutation inference
%   - Optionally tests an omnibus repeated-measures condition effect
%   - Computes map strength, weighted centroid, spatial dispersion, and
%     hemispheric lateralization for every participant and condition
%   - Tests the summary measures with paired tests and FDR correction
%   - Optionally generates spatial and summary figures
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - BROADNESS                         : Structure outputted by
%                                       BROADNESS_NetworkEstimation. It must
%                                       contain participant spatial patterns:
%      - .ActivationPatterns_BrainNetworks_Individual
%                                       sources x components x conditions
%                                       x participants
%
%  - Optional arguments (name-value pairs):
%      - 'principalcomps'              : Components to test (default: [1 2])
%      - 'conditionpairs'              : Rows containing condition pairs;
%                                       default: all condition pairs
%      - 'conditionnames'              : One label per condition
%      - 'mni_coords'                  : Source MNI coordinates (sources x 3).
%                                       The bundled 8-mm coordinates are
%                                       used when dimensions match
%      - 'mapnormalization'            : 'none' or 'rms' (default: 'none')
%      - 'omnibus'                     : 'on' or 'off' (default: 'off')
%      - 'permutations'                : Number of permutations (default: 1000)
%      - 'clusteralpha'                : Cluster-forming alpha (default: 0.01)
%      - 'alpha'                       : Corrected alpha (default: 0.05)
%      - 'neighbourdistance'           : Maximum MNI distance between spatial
%                                       neighbours. Default: inferred grid step
%      - 'randomseed'                  : Reproducible random seed (default: 1)
%      - 'figuremode'                  : 'off', 'show', 'save', or 'both'
%                                       (default: 'off')
%      - 'figurelayout'                : 'individual', 'summary', or 'both'
%                                       (default: 'both')
%      - 'OutputPath'                  : Base folder for saved figures
%      - 'figureformats'               : 'png', 'pdf', 'fig', or a cell array
%      - 'figureprefix'                : Optional saved-file prefix
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - SPATIAL_STATS                     : Structure containing:
%      - .Voxelwise.Pairwise           : Difference, t, effect-size, cluster,
%                                       corrected-p, and significance maps
%      - .Voxelwise.Omnibus            : Optional F and cluster results
%      - .Summary                      : Participant-level spatial measures
%      - .SummaryTests                 : Paired tests and omnibus ANOVA
%      - .Settings                     : Complete analysis settings
%  - FIGURES                           : Figure handles and saved paths
%
% ------------------------------------------------------------------------
%  NOTES:
% ------------------------------------------------------------------------
%  - Tests are paired: participants must occur in the same order in every
%    condition.
%  - Cluster correction is two-sided for paired contrasts and controls the
%    family-wise error rate across sources and selected components within
%    each condition contrast. Contrasts are treated as separate planned
%    statistical families.
%  - A significant cluster supports inference at the cluster level, not at
%    each individual source inside the cluster.
%  - The cluster-forming threshold is applied to the statistical map, not
%    to the original spatial activation patterns.
%  - Summary-measure FDR correction is applied jointly across measures,
%    selected components, and condition pairs.
%  - 'rms' normalization tests spatial redistribution after removing each
%    map's overall RMS strength. Summary measures are always calculated
%    from the original, unnormalized maps.
%  - The lateralization index is (left-right)/(left+right), calculated from
%    absolute weights. Sources at x = 0 are excluded.
%  - The default 1000 permutations is intended for practical exploration.
%    Use at least 5000, and preferably 10000, for final inference.
%
% ========================================================================
%  AUTHORS:
%  Leonardo Bonetti, Mathias Houe Andersen, Mattia Rosso
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  mathias.houe.andersen@regionh.dk
%  mattia.rosso@clin.au.dk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Danish Research Centre for Magnetic Resonance, Copenhagen University Hospital
%  Faculty of Health and Medical Sciences, University of Copenhagen
%  Aarhus (DK), Copenhagen (DK), Oxford (UK)
%
% ========================================================================

%% ----------------------------- Parse inputs -----------------------------

disp('Checking spatial-pattern statistical inputs')

opts = struct( ...
    'principalcomps', [1 2], ...
    'conditionpairs', [], ...
    'conditionnames', {{}}, ...
    'mni_coords', [], ...
    'mapnormalization', 'none', ...
    'omnibus', 'off', ...
    'permutations', 1000, ...
    'clusteralpha', 0.01, ...
    'alpha', 0.05, ...
    'neighbourdistance', [], ...
    'randomseed', 1, ...
    'figuremode', 'off', ...
    'figurelayout', 'both', ...
    'outputpath', [], ...
    'figureformats', {{'png'}}, ...
    'figureprefix', '');
opts = parse_name_value_pairs(opts, varargin{:});

if ~isstruct(BROADNESS) || ...
        ~isfield(BROADNESS,'ActivationPatterns_BrainNetworks_Individual')
    error(['BROADNESS must contain ' ...
        '"ActivationPatterns_BrainNetworks_Individual".']);
end

individualPatterns = BROADNESS.ActivationPatterns_BrainNetworks_Individual;
if ~isnumeric(individualPatterns) || isempty(individualPatterns) || ...
        ndims(individualPatterns) ~= 4
    error(['Participant spatial patterns must be a non-empty numeric array ' ...
        'with dimensions sources x components x conditions x participants.']);
end

[numberSources,numberComponents,numberConditions,numberParticipants] = ...
    size(individualPatterns);
if numberConditions < 2 || numberParticipants < 2
    error('At least two conditions and two participants are required.');
end

selectedComponents = opts.principalcomps(:)';
if ~isnumeric(selectedComponents) || isempty(selectedComponents) || ...
        any(~isfinite(selectedComponents)) || any(selectedComponents < 1) || ...
        any(selectedComponents > numberComponents) || ...
        any(fix(selectedComponents) ~= selectedComponents) || ...
        length(unique(selectedComponents)) ~= length(selectedComponents)
    error('"principalcomps" must contain unique valid component indices.');
end

if isempty(opts.conditionpairs)
    conditionPairs = nchoosek(1:numberConditions,2);
else
    conditionPairs = opts.conditionpairs;
    if ~isnumeric(conditionPairs) || size(conditionPairs,2) ~= 2 || ...
            isempty(conditionPairs) || any(~isfinite(conditionPairs(:))) || ...
            any(conditionPairs(:) < 1) || ...
            any(conditionPairs(:) > numberConditions) || ...
            any(fix(conditionPairs(:)) ~= conditionPairs(:)) || ...
            any(conditionPairs(:,1) == conditionPairs(:,2))
        error('"conditionpairs" must contain valid pairs of different conditions.');
    end
end

conditionNames = normalize_names(opts.conditionnames,numberConditions,'Condition');
mapNormalization = lower(char(string(opts.mapnormalization)));
if ~ismember(mapNormalization,{'none','rms'})
    error('"mapnormalization" must be ''none'' or ''rms''.');
end
validate_on_off(opts.omnibus,'omnibus');
validate_probability(opts.alpha,'alpha');
validate_probability(opts.clusteralpha,'clusteralpha');
if opts.clusteralpha >= opts.alpha
    error('"clusteralpha" must be smaller than "alpha".');
end
if ~isnumeric(opts.permutations) || ~isscalar(opts.permutations) || ...
        ~isfinite(opts.permutations) || opts.permutations < 3 || ...
        fix(opts.permutations) ~= opts.permutations
    error('"permutations" must be an integer of at least 3.');
end
if ~isnumeric(opts.randomseed) || ~isscalar(opts.randomseed) || ...
        ~isfinite(opts.randomseed) || fix(opts.randomseed) ~= opts.randomseed
    error('"randomseed" must be a finite integer scalar.');
end

coordinates = obtain_coordinates(opts.mni_coords,numberSources);
validSources = all(isfinite(coordinates),2);
if sum(validSources) < 2
    error('At least two sources must have valid MNI coordinates.');
end
if any(~isfinite(individualPatterns(validSources,selectedComponents,:,:)),'all')
    error('Spatial activation patterns must be finite at all valid MNI sources.');
end

if isempty(opts.neighbourdistance)
    neighbourDistance = infer_grid_step(coordinates(validSources,:));
else
    neighbourDistance = opts.neighbourdistance;
    if ~isnumeric(neighbourDistance) || ~isscalar(neighbourDistance) || ...
            ~isfinite(neighbourDistance) || neighbourDistance <= 0
        error('"neighbourdistance" must be a positive numeric scalar.');
    end
end

neighbours = build_neighbours(coordinates(validSources,:),neighbourDistance);
if ~any(cellfun(@(x) ~isempty(x),neighbours))
    error(['No spatial neighbours were found. Check "mni_coords" or ' ...
        'increase "neighbourdistance".']);
end

figureSettings = BROADNESS_FigureSettings(opts.figuremode,opts.figurelayout, ...
    opts.outputpath,opts.figureformats,opts.figureprefix, ...
    'SpatialPatternStatistics');
figureHandles = gobjects(0);
figureFiles = {};

%% -------------------------- Prepare spatial maps ------------------------

disp('Preparing unthresholded participant spatial activation patterns')

rawMaps = individualPatterns(:,selectedComponents,:,:);
analysisMaps = rawMaps;
if strcmp(mapNormalization,'rms')
    for component = 1:length(selectedComponents)
        for condition = 1:numberConditions
            for participant = 1:numberParticipants
                currentMap = rawMaps(validSources,component,condition,participant);
                mapRMS = sqrt(mean(currentMap.^2));
                if mapRMS > 0
                    analysisMaps(:,component,condition,participant) = ...
                        rawMaps(:,component,condition,participant) ./ mapRMS;
                else
                    error(['RMS normalization cannot be applied to an ' ...
                        'all-zero spatial map.']);
                end
            end
        end
    end
end

rngState = rng;
restoreRng = onCleanup(@() rng(rngState));
rng(opts.randomseed,'twister');

%% ------------------- Paired voxelwise cluster statistics ----------------

disp('Computing paired spatial cluster-permutation tests')

numberSelectedComponents = length(selectedComponents);
numberPairs = size(conditionPairs,1);
meanDifference = nan(numberSources,numberSelectedComponents,numberPairs);
tStatistics = nan(size(meanDifference));
cohenDz = nan(size(meanDifference));
uncorrectedP = nan(size(meanDifference));
clusterLabels = zeros(size(meanDifference));
clusterPMap = nan(size(meanDifference));
significant = false(size(meanDifference));
nullMaximumMass = nan(opts.permutations,numberPairs);
clusterResults = cell(numberSelectedComponents,numberPairs);

degreesFreedom = numberParticipants-1;
tThreshold = t_inverse_probability(1-opts.clusteralpha/2,degreesFreedom);

for pair = 1:numberPairs
    condition1 = conditionPairs(pair,1);
    condition2 = conditionPairs(pair,2);
    randomSigns = 2*(rand(numberParticipants,opts.permutations) > 0.5)-1;
    maximumMass = zeros(opts.permutations,1);

    differenceByComponent = cell(numberSelectedComponents,1);
    for component = 1:numberSelectedComponents
        differences = reshape(analysisMaps(validSources,component,condition1,:) - ...
            analysisMaps(validSources,component,condition2,:), ...
            sum(validSources),numberParticipants);
        differenceByComponent{component} = differences;

        [observedT,observedMean,observedDz] = paired_t_map(differences);
        sourceIndices = find(validSources);
        meanDifference(sourceIndices,component,pair) = observedMean;
        tStatistics(sourceIndices,component,pair) = observedT;
        cohenDz(sourceIndices,component,pair) = observedDz;
        uncorrectedP(sourceIndices,component,pair) = ...
            t_two_sided_probability(observedT,degreesFreedom);
    end

    % The matrix calculation makes sign-flip permutations fast. Spatial
    % clusters are then found separately in each permuted statistic map.
    for component = 1:numberSelectedComponents
        differences = differenceByComponent{component};
        sumSquares = sum(differences.^2,2);
        permutedMeans = (differences*randomSigns) ./ numberParticipants;
        permutedVariances = bsxfun(@minus,sumSquares, ...
            numberParticipants*permutedMeans.^2) ./ degreesFreedom;
        permutedT = permutedMeans ./ sqrt(max(permutedVariances,0) ./ ...
            numberParticipants);
        permutedT(~isfinite(permutedT)) = 0;

        for permutation = 1:opts.permutations
            [~,masses] = label_statistical_clusters( ...
                permutedT(:,permutation),tThreshold,neighbours,true);
            if ~isempty(masses)
                maximumMass(permutation) = max(maximumMass(permutation),max(masses));
            end
        end
    end
    nullMaximumMass(:,pair) = maximumMass;

    for component = 1:numberSelectedComponents
        observedT = tStatistics(validSources,component,pair);
        [labels,masses,clusterSigns] = label_statistical_clusters( ...
            observedT,tThreshold,neighbours,true);
        sourceIndices = find(validSources);
        clusterLabels(sourceIndices,component,pair) = labels;
        currentClusters = empty_cluster_structure();

        for cluster = 1:length(masses)
            correctedP = (1+sum(maximumMass >= masses(cluster))) / ...
                (opts.permutations+1);
            localMembers = find(labels == cluster);
            members = sourceIndices(localMembers);
            clusterPMap(members,component,pair) = correctedP;
            significant(members,component,pair) = correctedP <= opts.alpha;

            currentClusters(cluster).Cluster = cluster;
            currentClusters(cluster).Sign = clusterSigns(cluster);
            currentClusters(cluster).Mass = masses(cluster);
            currentClusters(cluster).CorrectedP = correctedP;
            currentClusters(cluster).Significant = correctedP <= opts.alpha;
            currentClusters(cluster).SourceIndices = members;
            currentClusters(cluster).MNI_Coordinates = coordinates(members,:);
        end
        clusterResults{component,pair} = currentClusters;
    end

    disp(['Condition contrast ' num2str(pair) ' / ' num2str(numberPairs) ...
        ' completed (' conditionNames{condition1} ' vs ' ...
        conditionNames{condition2} ')'])
end

contrastLabels = cell(numberPairs,1);
for pair = 1:numberPairs
    contrastLabels{pair} = [conditionNames{conditionPairs(pair,1)} ...
        ' vs ' conditionNames{conditionPairs(pair,2)}];
end

%% ---------------- Optional omnibus voxelwise condition test -------------

omnibus = struct();
omnibus.Enabled = strcmpi(opts.omnibus,'on');
if omnibus.Enabled
    disp('Computing omnibus spatial cluster-permutation test')
    dfCondition = numberConditions-1;
    dfError = (numberParticipants-1)*dfCondition;
    fThreshold = f_inverse_probability(1-opts.clusteralpha, ...
        dfCondition,dfError);
    fStatistics = nan(numberSources,numberSelectedComponents);
    omnibusP = nan(size(fStatistics));
    omnibusLabels = zeros(size(fStatistics));
    omnibusClusterPMap = nan(size(fStatistics));
    omnibusSignificant = false(size(fStatistics));
    omnibusClusters = cell(numberSelectedComponents,1);
    omnibusMaximumMass = zeros(opts.permutations,1);

    permutationOrders = zeros(numberConditions,numberParticipants,opts.permutations);
    for permutation = 1:opts.permutations
        for participant = 1:numberParticipants
            permutationOrders(:,participant,permutation) = randperm(numberConditions);
        end
    end

    for component = 1:numberSelectedComponents
        componentMaps = reshape(analysisMaps(validSources,component,:,:), ...
            sum(validSources),numberConditions,numberParticipants);
        observedF = repeated_measures_f_map(componentMaps);
        sourceIndices = find(validSources);
        fStatistics(sourceIndices,component) = observedF;
        omnibusP(sourceIndices,component) = ...
            f_upper_probability(observedF,dfCondition,dfError);

        for permutation = 1:opts.permutations
            permutedMaps = componentMaps;
            for participant = 1:numberParticipants
                permutedMaps(:,:,participant) = componentMaps(:, ...
                    permutationOrders(:,participant,permutation),participant);
            end
            permutedF = repeated_measures_f_map(permutedMaps);
            [~,masses] = label_statistical_clusters( ...
                permutedF,fThreshold,neighbours,false);
            if ~isempty(masses)
                omnibusMaximumMass(permutation) = max( ...
                    omnibusMaximumMass(permutation),max(masses));
            end
        end
    end

    for component = 1:numberSelectedComponents
        [labels,masses] = label_statistical_clusters( ...
            fStatistics(validSources,component),fThreshold,neighbours,false);
        sourceIndices = find(validSources);
        omnibusLabels(sourceIndices,component) = labels;
        currentClusters = empty_cluster_structure();
        for cluster = 1:length(masses)
            correctedP = (1+sum(omnibusMaximumMass >= masses(cluster))) / ...
                (opts.permutations+1);
            localMembers = find(labels == cluster);
            members = sourceIndices(localMembers);
            omnibusClusterPMap(members,component) = correctedP;
            omnibusSignificant(members,component) = correctedP <= opts.alpha;
            currentClusters(cluster).Cluster = cluster;
            currentClusters(cluster).Sign = 1;
            currentClusters(cluster).Mass = masses(cluster);
            currentClusters(cluster).CorrectedP = correctedP;
            currentClusters(cluster).Significant = correctedP <= opts.alpha;
            currentClusters(cluster).SourceIndices = members;
            currentClusters(cluster).MNI_Coordinates = coordinates(members,:);
        end
        omnibusClusters{component} = currentClusters;
    end

    omnibus.FStatistics = fStatistics;
    omnibus.PValues = omnibusP;
    omnibus.ClusterLabels = omnibusLabels;
    omnibus.ClusterCorrectedP = omnibusClusterPMap;
    omnibus.Significant = omnibusSignificant;
    omnibus.Clusters = omnibusClusters;
    omnibus.NullMaximumClusterMass = omnibusMaximumMass;
    omnibus.ClusterFormingThreshold = fThreshold;
    omnibus.DegreesFreedom = [dfCondition dfError];
    omnibus.CorrectionFamily = ...
        'All valid sources and selected components';
else
    omnibus.Reason = 'The optional omnibus condition test was switched off.';
end

%% -------------------- Participant spatial summary measures --------------

disp('Computing participant spatial summary measures')

summaryNames = {'MapStrength','CentroidX','CentroidY','CentroidZ', ...
    'SpatialDispersion','Lateralization'};
numberMeasures = length(summaryNames);
summaryValues = nan(numberMeasures,numberSelectedComponents, ...
    numberConditions,numberParticipants);

validCoordinates = coordinates(validSources,:);
leftMask = validCoordinates(:,1) < 0;
rightMask = validCoordinates(:,1) > 0;

for component = 1:numberSelectedComponents
    for condition = 1:numberConditions
        for participant = 1:numberParticipants
            spatialMap = rawMaps(validSources,component,condition,participant);
            absoluteWeights = abs(spatialMap);
            totalWeight = sum(absoluteWeights);
            mapStrength = sqrt(mean(spatialMap.^2));

            if totalWeight > 0
                centroid = sum(bsxfun(@times,validCoordinates,absoluteWeights),1) ./ ...
                    totalWeight;
                squaredDistance = sum(bsxfun(@minus,validCoordinates,centroid).^2,2);
                spatialDispersion = sqrt(sum(absoluteWeights.*squaredDistance) ./ ...
                    totalWeight);
            else
                centroid = [NaN NaN NaN];
                spatialDispersion = NaN;
            end

            leftWeight = sum(absoluteWeights(leftMask));
            rightWeight = sum(absoluteWeights(rightMask));
            if leftWeight+rightWeight > 0
                lateralization = (leftWeight-rightWeight) / ...
                    (leftWeight+rightWeight);
            else
                lateralization = NaN;
            end

            summaryValues(:,component,condition,participant) = ...
                [mapStrength centroid spatialDispersion lateralization];
        end
    end
end

summaryP = nan(numberMeasures,numberSelectedComponents,numberPairs);
summaryT = nan(size(summaryP));
summaryDz = nan(size(summaryP));
summaryMeanDifference = nan(size(summaryP));

for measure = 1:numberMeasures
    for component = 1:numberSelectedComponents
        for pair = 1:numberPairs
            values1 = squeeze(summaryValues(measure,component, ...
                conditionPairs(pair,1),:));
            values2 = squeeze(summaryValues(measure,component, ...
                conditionPairs(pair,2),:));
            differences = values1-values2;
            validDifferences = differences(isfinite(differences));
            numberValid = length(validDifferences);
            if numberValid < 2
                continue
            end
            differenceMean = mean(validDifferences);
            differenceSD = std(validDifferences,0);
            if differenceSD > 0
                currentT = differenceMean/(differenceSD/sqrt(numberValid));
                currentP = t_two_sided_probability(currentT,numberValid-1);
            elseif differenceMean ~= 0
                currentT = sign(differenceMean)*Inf;
                currentP = 0;
            else
                currentT = 0;
                currentP = 1;
            end
            summaryP(measure,component,pair) = currentP;
            summaryT(measure,component,pair) = currentT;
            summaryMeanDifference(measure,component,pair) = ...
                differenceMean;
            summaryDz(measure,component,pair) = ...
                differenceMean ./ differenceSD;
        end
    end
end
summaryDz(~isfinite(summaryDz)) = NaN;
[summarySignificant,summaryAdjustedP,summaryCriticalP] = ...
    BROADNESS_FDRCorrection(summaryP,opts.alpha);

if strcmpi(opts.omnibus,'on')
    summaryANOVA = BROADNESS_RepeatedMeasuresANOVA(summaryValues,opts.alpha);
else
    summaryANOVA.Enabled = false;
    summaryANOVA.Reason = 'The optional omnibus condition test was switched off.';
end

%% ------------------------------ Figures ---------------------------------

if ~strcmpi(opts.figuremode,'off')
    if figureSettings.MakeIndividual
        [spatialHandles,spatialFiles] = plot_spatial_results( ...
            coordinates,validSources,tStatistics,significant,selectedComponents, ...
            conditionPairs,contrastLabels,figureSettings);
        figureHandles = [figureHandles; spatialHandles];
        figureFiles = [figureFiles; spatialFiles];
    end
    if figureSettings.MakeSummary
        [summaryHandles,summaryFiles] = plot_summary_results( ...
            summaryValues,summaryNames,selectedComponents,conditionNames, ...
            figureSettings);
        figureHandles = [figureHandles; summaryHandles];
        figureFiles = [figureFiles; summaryFiles];
    end
end

FIGURES.Handles = figureHandles;
FIGURES.Files = figureFiles;

%% ---------------------------- Output structure --------------------------

SPATIAL_STATS = struct();
SPATIAL_STATS.MNI_Coordinates = coordinates;
SPATIAL_STATS.ValidSources = validSources;

SPATIAL_STATS.Voxelwise.Pairwise.ConditionPairs = conditionPairs;
SPATIAL_STATS.Voxelwise.Pairwise.ContrastLabels = contrastLabels;
SPATIAL_STATS.Voxelwise.Pairwise.MeanDifference = meanDifference;
SPATIAL_STATS.Voxelwise.Pairwise.TStatistics = tStatistics;
SPATIAL_STATS.Voxelwise.Pairwise.CohenDz = cohenDz;
SPATIAL_STATS.Voxelwise.Pairwise.PValues = uncorrectedP;
SPATIAL_STATS.Voxelwise.Pairwise.ClusterLabels = clusterLabels;
SPATIAL_STATS.Voxelwise.Pairwise.ClusterCorrectedP = clusterPMap;
SPATIAL_STATS.Voxelwise.Pairwise.Significant = significant;
SPATIAL_STATS.Voxelwise.Pairwise.Clusters = clusterResults;
SPATIAL_STATS.Voxelwise.Pairwise.NullMaximumClusterMass = nullMaximumMass;
SPATIAL_STATS.Voxelwise.Pairwise.ClusterFormingThreshold = tThreshold;
SPATIAL_STATS.Voxelwise.Pairwise.DegreesFreedom = degreesFreedom;
SPATIAL_STATS.Voxelwise.Pairwise.CorrectionFamily = ...
    'All valid sources and selected components, separately per contrast';
SPATIAL_STATS.Voxelwise.Omnibus = omnibus;

SPATIAL_STATS.Summary.Names = summaryNames;
SPATIAL_STATS.Summary.Values = summaryValues;
for measure = 1:numberMeasures
    SPATIAL_STATS.Summary.(summaryNames{measure}) = ...
        squeeze(summaryValues(measure,:,:,:));
end

SPATIAL_STATS.SummaryTests.Pairwise.ConditionPairs = conditionPairs;
SPATIAL_STATS.SummaryTests.Pairwise.ContrastLabels = contrastLabels;
SPATIAL_STATS.SummaryTests.Pairwise.MeanDifference = summaryMeanDifference;
SPATIAL_STATS.SummaryTests.Pairwise.TStatistics = summaryT;
SPATIAL_STATS.SummaryTests.Pairwise.CohenDz = summaryDz;
SPATIAL_STATS.SummaryTests.Pairwise.PValues = summaryP;
SPATIAL_STATS.SummaryTests.Pairwise.AdjustedPValues = summaryAdjustedP;
SPATIAL_STATS.SummaryTests.Pairwise.Significant = summarySignificant;
SPATIAL_STATS.SummaryTests.Pairwise.FDRCriticalP = summaryCriticalP;
SPATIAL_STATS.SummaryTests.Pairwise.FDRMethod = 'Benjamini-Hochberg';
SPATIAL_STATS.SummaryTests.Pairwise.FDRFamily = ...
    'All summary measures, selected components, and condition pairs';
SPATIAL_STATS.SummaryTests.ANOVA = summaryANOVA;

SPATIAL_STATS.Settings.SelectedComponents = selectedComponents;
SPATIAL_STATS.Settings.ConditionPairs = conditionPairs;
SPATIAL_STATS.Settings.ConditionNames = conditionNames;
SPATIAL_STATS.Settings.MapNormalization = mapNormalization;
SPATIAL_STATS.Settings.Permutations = opts.permutations;
SPATIAL_STATS.Settings.ClusterFormingAlpha = opts.clusteralpha;
SPATIAL_STATS.Settings.Alpha = opts.alpha;
SPATIAL_STATS.Settings.NeighbourDistance = neighbourDistance;
SPATIAL_STATS.Settings.RandomSeed = opts.randomseed;
SPATIAL_STATS.Settings.Omnibus = lower(char(opts.omnibus));
SPATIAL_STATS.Settings.InferenceData = 'Unthresholded individual spatial patterns';
SPATIAL_STATS.Settings.ClusterCorrection = 'Maximum cluster mass';
SPATIAL_STATS.Figures = FIGURES;

disp('Spatial activation pattern statistics completed')

end


%% ========================================================================
%                              Helper functions
% =========================================================================

function [tValues,meanValues,cohenDz] = paired_t_map(differences)
numberParticipants = size(differences,2);
meanValues = mean(differences,2);
standardDeviation = std(differences,0,2);
tValues = meanValues ./ (standardDeviation ./ sqrt(numberParticipants));
cohenDz = meanValues ./ standardDeviation;
tValues(~isfinite(tValues)) = 0;
cohenDz(~isfinite(cohenDz)) = NaN;
end


function fValues = repeated_measures_f_map(data)
% Data dimensions: sources x conditions x participants.
[~,numberConditions,numberParticipants] = size(data);
grandMean = mean(mean(data,3),2);
conditionMeans = mean(data,3);
participantMeans = mean(data,2);
conditionSS = numberParticipants * sum(bsxfun(@minus,conditionMeans,grandMean).^2,2);
participantSS = numberConditions * sum(bsxfun(@minus,participantMeans,grandMean).^2,3);
totalSS = sum(sum(bsxfun(@minus,data,grandMean).^2,3),2);
errorSS = max(0,totalSS-conditionSS-participantSS);
dfCondition = numberConditions-1;
dfError = (numberParticipants-1)*dfCondition;
fValues = (conditionSS./dfCondition) ./ (errorSS./dfError);
fValues(~isfinite(fValues)) = 0;
fValues = fValues(:);
end


function probability = t_two_sided_probability(tValue,degreesFreedom)
probability = nan(size(tValue));
finiteValues = isfinite(tValue);
betaValue = degreesFreedom ./ ...
    (degreesFreedom+tValue(finiteValues).^2);
probability(finiteValues) = betainc(betaValue,degreesFreedom/2,0.5);
probability(isinf(tValue)) = 0;
end


function tValue = t_inverse_probability(probability,degreesFreedom)
% Positive Student-t quantile for probability values greater than 0.5.
betaValue = betaincinv(2*(1-probability),degreesFreedom/2,0.5);
tValue = sqrt(degreesFreedom*(1-betaValue)/betaValue);
end


function probability = f_upper_probability(fValue,dfNumerator,dfDenominator)
probability = nan(size(fValue));
finiteValues = isfinite(fValue) & fValue >= 0;
betaValue = (dfNumerator*fValue(finiteValues)) ./ ...
    (dfNumerator*fValue(finiteValues)+dfDenominator);
probability(finiteValues) = betainc(betaValue,dfNumerator/2, ...
    dfDenominator/2,'upper');
probability(isinf(fValue) & fValue > 0) = 0;
end


function fValue = f_inverse_probability(probability,dfNumerator,dfDenominator)
betaValue = betaincinv(probability,dfNumerator/2,dfDenominator/2);
fValue = (dfDenominator*betaValue) / ...
    (dfNumerator*(1-betaValue));
end


function [labels,masses,clusterSigns] = label_statistical_clusters( ...
        statisticMap,threshold,neighbours,twoSided)
statisticMap = statisticMap(:);
labels = zeros(size(statisticMap));
masses = [];
clusterSigns = [];

if twoSided
    masks = {statisticMap >= threshold, statisticMap <= -threshold};
    signs = [1 -1];
else
    masks = {statisticMap >= threshold};
    signs = 1;
end

clusterNumber = 0;
for signIndex = 1:length(masks)
    active = masks{signIndex};
    visited = false(size(active));
    seeds = find(active);
    for seedIndex = 1:length(seeds)
        seed = seeds(seedIndex);
        if visited(seed)
            continue
        end
        clusterNumber = clusterNumber+1;
        queue = seed;
        visited(seed) = true;
        members = zeros(sum(active),1);
        memberCount = 0;
        queueStart = 1;
        while queueStart <= length(queue)
            current = queue(queueStart);
            queueStart = queueStart+1;
            memberCount = memberCount+1;
            members(memberCount) = current;
            currentNeighbours = neighbours{current};
            newNeighbours = currentNeighbours(active(currentNeighbours) & ...
                ~visited(currentNeighbours));
            if ~isempty(newNeighbours)
                visited(newNeighbours) = true;
                queue = [queue; newNeighbours(:)]; %#ok<AGROW>
            end
        end
        members = members(1:memberCount);
        labels(members) = clusterNumber;
        masses(clusterNumber,1) = sum(abs(statisticMap(members))); %#ok<AGROW>
        clusterSigns(clusterNumber,1) = signs(signIndex); %#ok<AGROW>
    end
end
end


function neighbours = build_neighbours(coordinates,neighbourDistance)
% Calculate distances in small blocks to avoid requiring an additional
% neighbour-search function or allocating one full sources x sources matrix.
numberSources = size(coordinates,1);
neighbours = cell(numberSources,1);
squaredNorms = sum(coordinates.^2,2)';
maximumSquaredDistance = (neighbourDistance*1.001)^2;
blockSize = 250;
for firstSource = 1:blockSize:numberSources
    lastSource = min(numberSources,firstSource+blockSize-1);
    blockCoordinates = coordinates(firstSource:lastSource,:);
    squaredDistances = bsxfun(@plus,sum(blockCoordinates.^2,2),squaredNorms) - ...
        2*(blockCoordinates*coordinates');
    for localSource = 1:size(blockCoordinates,1)
        source = firstSource+localSource-1;
        currentNeighbours = find(squaredDistances(localSource,:) <= ...
            maximumSquaredDistance & squaredDistances(localSource,:) > 1e-8);
        neighbours{source} = currentNeighbours(:);
    end
end
end


function gridStep = infer_grid_step(coordinates)
axisSteps = [];
for dimension = 1:3
    uniqueValues = unique(coordinates(:,dimension));
    differences = diff(uniqueValues);
    axisSteps = [axisSteps; differences(differences > 1e-6)]; %#ok<AGROW>
end
if isempty(axisSteps)
    error('The MNI grid spacing could not be inferred.');
end
gridStep = min(axisSteps);
end


function coordinates = obtain_coordinates(inputCoordinates,numberSources)
if isempty(inputCoordinates)
    coordinateFile = which('MNI152_8mm_coord_dyi.mat');
    if isempty(coordinateFile)
        error(['MNI coordinates were not provided and the bundled coordinate ' ...
            'file is not available on the MATLAB path.']);
    end
    coordinateStructure = load(coordinateFile);
    if ~isfield(coordinateStructure,'MNI8')
        error('The bundled coordinate file does not contain "MNI8".');
    end
    coordinates = coordinateStructure.MNI8;
else
    coordinates = inputCoordinates;
end
if ~isnumeric(coordinates) || size(coordinates,1) ~= numberSources || ...
        size(coordinates,2) ~= 3
    error('"mni_coords" must have dimensions sources x 3.');
end
end


function clusters = empty_cluster_structure()
clusters = struct('Cluster',{},'Sign',{},'Mass',{},'CorrectedP',{}, ...
    'Significant',{},'SourceIndices',{},'MNI_Coordinates',{});
end


function [handles,files] = plot_spatial_results(coordinates,validSources, ...
        tStatistics,significant,selectedComponents,conditionPairs, ...
        contrastLabels,figureSettings)
handles = gobjects(0);
files = {};
templateFigure = [];
templateAxes = [];
templateFile = which('BrainTemplate_MNI152_1mm_FullBrain.fig');
if ~isempty(templateFile)
    templateFigure = openfig(templateFile,'new','invisible');
    templateAxes = findobj(templateFigure,'Type','axes');
    if ~isempty(templateAxes)
        templateAxes = templateAxes(1);
    end
end

validCoordinates = coordinates(validSources,:);
for component = 1:length(selectedComponents)
    for pair = 1:size(conditionPairs,1)
        fig = figure('Visible',figureSettings.Visible,'Color','w', ...
            'Position',[100 100 1150 500]);
        layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

        ax1 = nexttile(layout);
        copy_brain_template(ax1,templateAxes)
        currentT = tStatistics(validSources,component,pair);
        maximumT = max(abs(currentT));
        if maximumT == 0 || ~isfinite(maximumT), maximumT = 1; end
        scatter3(ax1,validCoordinates(:,1),validCoordinates(:,2), ...
            validCoordinates(:,3),12,currentT,'filled', ...
            'MarkerFaceAlpha',0.65,'MarkerEdgeAlpha',0.25);
        colormap(ax1,blue_white_red(256))
        set(ax1,'CLim',[-maximumT maximumT])
        colorbar(ax1)
        title(ax1,'Paired t-statistic')
        finish_brain_axes(ax1)

        ax2 = nexttile(layout);
        copy_brain_template(ax2,templateAxes)
        currentSignificant = significant(validSources,component,pair);
        positive = currentSignificant & currentT > 0;
        negative = currentSignificant & currentT < 0;
        if any(positive)
            scatter3(ax2,validCoordinates(positive,1),validCoordinates(positive,2), ...
                validCoordinates(positive,3),22,[0.76 0.16 0.14],'filled');
        end
        if any(negative)
            scatter3(ax2,validCoordinates(negative,1),validCoordinates(negative,2), ...
                validCoordinates(negative,3),22,[0.14 0.32 0.72],'filled');
        end
        if ~any(currentSignificant)
            text(ax2,0.5,0.04,'No cluster survived correction', ...
                'Units','normalized','HorizontalAlignment','center', ...
                'FontAngle','italic','Color',[0.35 0.35 0.35]);
        end
        title(ax2,'Cluster-corrected differences')
        finish_brain_axes(ax2)

        title(layout,['Brain Network ' num2str(selectedComponents(component)) ...
            ': ' contrastLabels{pair}],'FontWeight','normal')
        rotate3d(fig,'on')
        handles(end+1,1) = fig; %#ok<AGROW>
        savedFiles = BROADNESS_FinalizeFigure(fig,figureSettings, ...
            ['SpatialStatistics_BN' num2str(selectedComponents(component)) ...
            '_Condition' num2str(conditionPairs(pair,1)) '_vs_' ...
            num2str(conditionPairs(pair,2))]);
        files = [files; savedFiles]; %#ok<AGROW>
    end
end
if ~isempty(templateFigure) && isgraphics(templateFigure,'figure')
    close(templateFigure)
end
end


function [handles,files] = plot_summary_results(summaryValues,summaryNames, ...
        selectedComponents,conditionNames,figureSettings)
handles = gobjects(0);
files = {};
numberConditions = length(conditionNames);
colors = lines(numberConditions);
yLabels = {'RMS weight','MNI x (mm)','MNI y (mm)','MNI z (mm)', ...
    'Distance (mm)','(Left - right) / (left + right)'};

for component = 1:length(selectedComponents)
    fig = figure('Visible',figureSettings.Visible,'Color','w', ...
        'Position',[100 80 1150 700]);
    layout = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    for measure = 1:length(summaryNames)
        ax = nexttile(layout);
        hold(ax,'on')
        values = squeeze(summaryValues(measure,component,:,:))';
        plot(ax,1:numberConditions,values','-','Color',[0.78 0.78 0.78], ...
            'LineWidth',0.5,'HandleVisibility','off');
        means = mean(values,1,'omitnan');
        standardErrors = std(values,0,1,'omitnan') ./ ...
            sqrt(sum(isfinite(values),1));
        for condition = 1:numberConditions
            scatter(ax,repmat(condition,size(values,1),1),values(:,condition), ...
                13,colors(condition,:),'filled','MarkerFaceAlpha',0.45, ...
                'HandleVisibility','off');
        end
        errorbar(ax,1:numberConditions,means,standardErrors,'k-', ...
            'LineWidth',1.4,'Marker','o','MarkerFaceColor','w');
        xlim(ax,[0.6 numberConditions+0.4])
        xticks(ax,1:numberConditions)
        xticklabels(ax,conditionNames)
        xtickangle(ax,25)
        ylabel(ax,yLabels{measure})
        title(ax,summaryNames{measure},'Interpreter','none')
        grid(ax,'on')
        box(ax,'off')
        set(ax,'FontSize',10,'LineWidth',1,'TickDir','out')
    end
    title(layout,['Spatial summaries - Brain Network ' ...
        num2str(selectedComponents(component))],'FontWeight','normal')
    handles(end+1,1) = fig; %#ok<AGROW>
    savedFiles = BROADNESS_FinalizeFigure(fig,figureSettings, ...
        ['SpatialSummaries_BN' num2str(selectedComponents(component))]);
    files = [files; savedFiles]; %#ok<AGROW>
end
end


function copy_brain_template(ax,templateAxes)
hold(ax,'on')
if ~isempty(templateAxes) && isgraphics(templateAxes,'axes')
    copyobj(allchild(templateAxes),ax);
    set(ax,'XLim',templateAxes.XLim,'YLim',templateAxes.YLim, ...
        'ZLim',templateAxes.ZLim,'View',templateAxes.View, ...
        'Projection',templateAxes.Projection, ...
        'DataAspectRatio',templateAxes.DataAspectRatio);
end
end


function finish_brain_axes(ax)
axis(ax,'equal')
axis(ax,'vis3d')
axis(ax,'off')
view(ax,[-90 10])
end


function colors = blue_white_red(numberColors)
half = floor(numberColors/2);
blue = [linspace(0.10,1,half)' linspace(0.25,1,half)' ones(half,1)];
redCount = numberColors-half;
red = [ones(redCount,1) linspace(1,0.18,redCount)' ...
    linspace(1,0.12,redCount)'];
colors = [blue; red];
end


function names = normalize_names(inputNames,numberNames,prefix)
if isempty(inputNames)
    names = arrayfun(@(x) [prefix ' ' num2str(x)],1:numberNames, ...
        'UniformOutput',false);
elseif isstring(inputNames)
    names = cellstr(inputNames(:));
elseif iscell(inputNames)
    names = inputNames(:);
else
    error('"conditionnames" must be a cell array or string array.');
end
if length(names) ~= numberNames
    error('Provide one entry in "conditionnames" for every condition.');
end
for index = 1:length(names)
    names{index} = char(string(names{index}));
end
end


function validate_on_off(value,name)
if ~(ischar(value) || (isstring(value) && isscalar(value))) || ...
        ~ismember(lower(char(value)),{'on','off'})
    error('"%s" must be ''on'' or ''off''.',name);
end
end


function validate_probability(value,name)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value <= 0 || value >= 1
    error('"%s" must be a numeric scalar between 0 and 1.',name);
end
end


function opts = parse_name_value_pairs(opts,varargin)
if mod(length(varargin),2) ~= 0
    error('Arguments must be given as name-value pairs.');
end
for index = 1:2:length(varargin)
    name = lower(char(string(varargin{index})));
    if ~isfield(opts,name)
        error('Unknown option "%s".',varargin{index});
    end
    opts.(name) = varargin{index+1};
end
end
