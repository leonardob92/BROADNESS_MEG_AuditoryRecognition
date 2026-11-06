function [significant, adjustedP, criticalP] = BROADNESS_FDRCorrection(pValues, alpha)

% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  FALSE DISCOVERY RATE CORRECTION
% ========================================================================
%
%  Applies the Benjamini-Hochberg false discovery rate (FDR) procedure to
%  the complete family of p-values provided by the user. NaN values are
%  retained as NaN and are excluded from the correction.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - pValues    : Numeric array containing uncorrected p-values
%  - alpha      : FDR level (default: 0.05)
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - significant: Logical array with the same dimensions as pValues
%  - adjustedP  : Benjamini-Hochberg adjusted p-values
%  - criticalP  : Largest uncorrected p-value passing the FDR criterion;
%                 NaN when no result passes
%
%  The user determines the correction family by choosing which p-values
%  are supplied together. The statistical examples pass all tested time
%  points, networks, and planned contrasts together.
%
% ========================================================================

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
if ~isnumeric(pValues) || ~isreal(pValues)
    error('"pValues" must be a real numeric array.');
end
if ~isnumeric(alpha) || ~isscalar(alpha) || ~isfinite(alpha) || ...
        alpha <= 0 || alpha >= 1
    error('"alpha" must be a numeric scalar between 0 and 1.');
end
if any(pValues(~isnan(pValues)) < 0 | pValues(~isnan(pValues)) > 1)
    error('Finite p-values must be between 0 and 1.');
end

originalSize = size(pValues);
pVector = pValues(:);
validIndices = find(~isnan(pVector));
adjustedVector = nan(size(pVector));
significantVector = false(size(pVector));
criticalP = NaN;

if isempty(validIndices)
    significant = reshape(significantVector, originalSize);
    adjustedP = reshape(adjustedVector, originalSize);
    return
end

[sortedP, sortingOrder] = sort(pVector(validIndices));
numberTests = length(sortedP);
fdrLimits = ((1:numberTests)' ./ numberTests) * alpha;
lastPassing = find(sortedP <= fdrLimits, 1, 'last');
if ~isempty(lastPassing)
    criticalP = sortedP(lastPassing);
    significantVector(validIndices) = pVector(validIndices) <= criticalP;
end

adjustedSorted = sortedP .* numberTests ./ (1:numberTests)';
for index = numberTests-1:-1:1
    adjustedSorted(index) = min(adjustedSorted(index), adjustedSorted(index+1));
end
adjustedSorted = min(adjustedSorted, 1);
inverseOrder = zeros(numberTests,1);
inverseOrder(sortingOrder) = 1:numberTests;
adjustedVector(validIndices) = adjustedSorted(inverseOrder);

significant = reshape(significantVector, originalSize);
adjustedP = reshape(adjustedVector, originalSize);

end
