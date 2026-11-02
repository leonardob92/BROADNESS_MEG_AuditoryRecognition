function ANOVA = BROADNESS_RepeatedMeasuresANOVA(data, alpha)

% ========================================================================
%  BROADNESS TIME-RESOLVED REPEATED-MEASURES ANOVA
% ========================================================================
%
%  Performs a one-factor repeated-measures ANOVA at every time-point and
%  network. Input data must have dimensions:
%  time x networks x conditions x participants.
%
%  For more than two conditions, Greenhouse-Geisser correction is applied
%  to the degrees of freedom. FDR is then applied jointly across all tested
%  time-points and networks.
%
% ========================================================================

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
if ~isnumeric(data) || ~isreal(data) || ndims(data) ~= 4
    error('"data" must be a real time x networks x conditions x participants array.');
end

[numberTimepoints,numberNetworks,numberConditions,~] = size(data);
if numberConditions < 2
    error('Repeated-measures ANOVA requires at least two conditions.');
end

fStatistics = nan(numberTimepoints,numberNetworks);
pValues = nan(size(fStatistics));
uncorrectedPValues = nan(size(fStatistics));
epsilonGG = nan(size(fStatistics));
conditionDF = nan(size(fStatistics));
errorDF = nan(size(fStatistics));

for timepoint = 1:numberTimepoints
    for network = 1:numberNetworks
        conditionData = reshape(data(timepoint,network,:,:), ...
            numberConditions,[])';
        conditionData = conditionData(all(isfinite(conditionData),2),:);
        numberParticipants = size(conditionData,1);
        if numberParticipants < 2
            continue
        end

        grandMean = mean(conditionData(:));
        conditionMeans = mean(conditionData,1);
        participantMeans = mean(conditionData,2);
        conditionSumSquares = numberParticipants * ...
            sum((conditionMeans-grandMean).^2);
        participantSumSquares = numberConditions * ...
            sum((participantMeans-grandMean).^2);
        totalSumSquares = sum((conditionData-grandMean).^2,'all');
        errorSumSquares = max(0,totalSumSquares-conditionSumSquares- ...
            participantSumSquares);

        dfCondition = numberConditions-1;
        dfError = (numberParticipants-1)*dfCondition;
        conditionMeanSquare = conditionSumSquares/dfCondition;
        errorMeanSquare = errorSumSquares/dfError;
        if errorMeanSquare == 0
            if conditionMeanSquare > 0
                fValue = Inf;
            else
                fValue = NaN;
            end
        else
            fValue = conditionMeanSquare/errorMeanSquare;
        end

        epsilon = 1;
        if numberConditions > 2
            centeringMatrix = eye(numberConditions) - ...
                ones(numberConditions)/numberConditions;
            covarianceMatrix = cov(conditionData);
            centeredCovariance = centeringMatrix*covarianceMatrix*centeringMatrix;
            epsilonDenominator = dfCondition * ...
                sum(centeredCovariance.^2,'all');
            if epsilonDenominator > 0
                epsilon = trace(centeredCovariance)^2/epsilonDenominator;
                epsilon = min(1,max(1/dfCondition,epsilon));
            end
        end

        fStatistics(timepoint,network) = fValue;
        epsilonGG(timepoint,network) = epsilon;
        conditionDF(timepoint,network) = epsilon*dfCondition;
        errorDF(timepoint,network) = epsilon*dfError;
        uncorrectedPValues(timepoint,network) = ...
            f_upper_probability(fValue,dfCondition,dfError);
        pValues(timepoint,network) = f_upper_probability( ...
            fValue,epsilon*dfCondition,epsilon*dfError);
    end
end

[significant,adjustedP,criticalP] = ...
    BROADNESS_FDRCorrection(pValues,alpha);

ANOVA.Enabled = true;
ANOVA.FStatistics = fStatistics;
ANOVA.PValues = pValues;
ANOVA.UncorrectedPValues = uncorrectedPValues;
ANOVA.AdjustedPValues = adjustedP;
ANOVA.Significant = significant;
ANOVA.FDRCriticalP = criticalP;
ANOVA.Alpha = alpha;
ANOVA.FDRMethod = 'Benjamini-Hochberg';
ANOVA.FDRFamily = 'All tested time-points and networks';
ANOVA.GreenhouseGeisserEpsilon = epsilonGG;
ANOVA.ConditionDF = conditionDF;
ANOVA.ErrorDF = errorDF;
ANOVA.Correction = 'Greenhouse-Geisser for more than two conditions';

end


function probability = f_upper_probability(fValue,dfNumerator,dfDenominator)

if isnan(fValue)
    probability = NaN;
elseif isinf(fValue)
    probability = 0;
else
    betaValue = (dfNumerator*fValue) / ...
        (dfNumerator*fValue+dfDenominator);
    probability = betainc(betaValue,dfNumerator/2,dfDenominator/2,'upper');
end

end
