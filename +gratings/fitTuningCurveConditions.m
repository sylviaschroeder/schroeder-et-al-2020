function [parameters, predictions] = ...
    fitTuningCurveConditions(responses, stimDirections, blanks, ...
    conditions, fixedPars, runFitTwice)

% responses         [stimuli x repetitions]; response strength for each
%                   trial
% stimDirections    [stimuli x 1]; direction
% blanks            [stimuli x 1]; true if blank, false otherwise
% conditions        [stimuli x repetitions]; integer for each trial,
%                   differentiation between different conditions, such as 
%                   running and not running
% fixedPars         [1 x p]; parameters that will be fixed across
%                   conditions, default: [1 5] -> pref. dir, width
% runFitTwice       if 0, all parameters are fit within given limits, if 1,
%                   parameters not fixed across conditions are fit a 2nd 
%                   time and the other parameters are set to values found
%                   in first run

% parameters        [parameters x conditions]; different orientation tuning
%                   parameters for each condition; parameters: (1)
%                   preferred direction, (2) ampl. at pref. direction, (3)
%                   direction selectivity index, (4) offset of tuning curve
%                   from zero, (5) tuning width
% predictions       [stimuli x repetitions]; model predictions

numFitIterations = 10;

if nargin < 6
    runFitTwice = 0;
end
if nargin < 5
    fixedPars = [1 5];
end

dirs = repmat(stimDirections,1,size(responses,2));
resp = responses(~blanks,:);
conds = conditions(~blanks,:);

ind = ~isnan(resp) & ~isnan(conds);
paramDeltas = NaN(1,5);
paramDeltas(fixedPars) = 0;
paramLimits = repmat([-Inf;Inf],1,5);
paramLimits(1,3) = 0; % limit min. DI to 0 (not -1) so that pref. dir. does not change
paramLimits(1,5) = median(diff(unique(dirs))); % limit min. sigma to the sample interval between tested orientations/directions
pars = gratings.fitoriConditions(dirs(ind), resp(ind), conds(ind), ...
    paramDeltas, paramLimits, numFitIterations);
if runFitTwice == 1 && ~isempty(fixedPars)
    paramLimits(:, fixedPars) = repmat(pars(fixedPars), 2, 1);
    pars = gratings.fitoriConditions(dirs(ind), resp(ind), conds(ind), ...
        [], paramLimits, numFitIterations);
end

parameters = reshape(pars',5,[]);
parameters(:,2:end) = bsxfun(@plus, parameters(:,1), parameters(:,2:end));

pred = NaN(size(resp));
pred(ind) = gratings.orituneWrappedConditions(pars, dirs(ind), conds(ind));

predictions = NaN(size(responses));
predictions(~blanks,:) = pred;