function error = fitTuningCurveConditions_forCrossVal(respTraining, ...
    respTest, variables, testIndices)
% respTraining  [trial x stimulus]; training set
% respTest      [trial x stimulus]; test set
% variables     {1 x 4}, {1}: [stimulus x 1], direction
%                        {2}: [stimulus x 1], true if blank, false
%                        otherwise
%                        {3}: [trial x stimulus], condition of each trial
%                        {4}: [1 x p], parameters that are fixed across
%                        conditions
% (testIndices) [n x 1], stimulus indices that are part of test set

% error         [n x 1], difference between prediction and data for each
%               entry in test set

directions = variables{1};
blanks = variables{2};
conditions = variables{3};
fixedPars = variables{4};

parameters = gratings.fitTuningCurveConditions(respTraining', directions, blanks, ...
    conditions', fixedPars, 1);

if nargin < 4
    testIndices = find(isnan(respTraining));
end
[~,stim] = ind2sub(size(respTraining), testIndices);

conds = conditions(testIndices);
dirs = NaN(length(testIndices),1);
for k = 1:length(testIndices)
    j = stim(k);
    if ~isempty(j)
        dirs(k) = directions(j);
    end
end
invalid = isnan(dirs);
conds(invalid) = [];
dirs(invalid) = [];

cs = unique(conds);
prediction = NaN(length(dirs),1);
for c = 1:length(cs)
    j = conds==cs(c);
    prediction(j) = gratings.orituneWrappedConditions(parameters(:,c), dirs(j), conds(j));
end
test = respTest(testIndices);
test(invalid) = [];
error = NaN(length(testIndices), 1);
error(~invalid) = prediction - test;