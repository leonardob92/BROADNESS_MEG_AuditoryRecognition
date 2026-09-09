function windows = BROADNESS_SignificantSamplesToWindows(significantSamples, time, minimumSamples)

% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  SIGNIFICANT SAMPLES TO TIME WINDOWS
% ========================================================================
%
%  Converts consecutive significant samples into the nested time-window
%  format accepted by BROADNESS_Plot_ActivationTimeseries.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - significantSamples : Logical or binary vector, one value per time-point
%  - time               : Corresponding time vector in seconds
%  - minimumSamples     : Minimum consecutive samples retained (default: 1)
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENT:
% ------------------------------------------------------------------------
%  - windows            : Cell array containing [start_time end_time] pairs
%
%  This utility only reformats an already corrected statistical result. It
%  does not perform an additional statistical or cluster-level correction.
%
% ========================================================================

if nargin < 3 || isempty(minimumSamples)
    minimumSamples = 1;
end
if ~isvector(significantSamples) || ~isvector(time) || ...
        numel(significantSamples) ~= numel(time)
    error('"significantSamples" and "time" must be vectors of equal length.');
end
if any(~isfinite(time)) || any(diff(time(:)) <= 0)
    error('"time" must contain finite, strictly increasing values.');
end
if ~isnumeric(minimumSamples) || ~isscalar(minimumSamples) || ...
        minimumSamples < 1 || fix(minimumSamples) ~= minimumSamples
    error('"minimumSamples" must be a positive integer.');
end

significantSamples = logical(significantSamples(:));
changes = diff([false; significantSamples; false]);
starts = find(changes == 1);
ends = find(changes == -1) - 1;
keep = (ends - starts + 1) >= minimumSamples;
starts = starts(keep);
ends = ends(keep);

windows = cell(length(starts),1);
for window = 1:length(starts)
    windows{window} = [time(starts(window)) time(ends(window))];
end

end
