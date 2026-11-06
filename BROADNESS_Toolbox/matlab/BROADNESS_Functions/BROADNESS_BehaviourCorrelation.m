function [correlations, pValues, adjustedP, significant, criticalP] = ...
    BROADNESS_BehaviourCorrelation(data, behaviour, alpha, correlationType)

% ========================================================================
%  BROADNESS CORRELATION WITH BEHAVIOUR
% ========================================================================
%
%  Correlates every measure in "data" with one behavioural value per
%  participant and applies Benjamini-Hochberg FDR correction. Participants
%  must be stored in the final dimension of "data".
%
%  INPUTS:
%  - data            : Numeric array with participants in the last dimension
%  - behaviour       : One value per participant
%  - alpha           : FDR level (default: 0.05)
%  - correlationType : 'Pearson' (default) or 'Spearman'
%
%  OUTPUTS have the same dimensions as "data", excluding participants.
%
% ========================================================================

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 4 || isempty(correlationType)
    correlationType = 'Pearson';
end
if ~isnumeric(data) || ~isreal(data) || ndims(data) < 2
    error('"data" must be a real numeric array with participants last.');
end
if ~isnumeric(behaviour) || ~isvector(behaviour) || ...
        numel(behaviour) ~= size(data,ndims(data))
    error('"behaviour" must contain one numeric value per participant.');
end
if ~ismember(lower(correlationType), {'pearson','spearman'})
    error('"correlationType" must be ''Pearson'' or ''Spearman''.');
end

behaviour = behaviour(:);
if sum(isfinite(behaviour)) < 3
    error('"behaviour" must contain at least three finite participant values.');
end

dataSize = size(data);
numberParticipants = dataSize(end);
outputSize = dataSize(1:end-1);
if isscalar(outputSize)
    outputSize = [outputSize 1];
end
participantData = reshape(data,[],numberParticipants)';

[correlationVector,pValueVector] = corr(participantData,behaviour, ...
    'Type',correlationType,'Rows','pairwise');

correlations = reshape(correlationVector,outputSize);
pValues = reshape(pValueVector,outputSize);
[significant,adjustedP,criticalP] = BROADNESS_FDRCorrection(pValues,alpha);

end
