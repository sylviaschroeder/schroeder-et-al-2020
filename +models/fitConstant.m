function error = fitConstant(respTraining, respTest, variables, testIndices)
% respTraining  [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      {1 x 3}, {1}: [stimulus x 2], 1st col: direction, 2nd col:
%                             stimID
%                        {2}: [n x 1], stim IDs of blank stimuli, n is
%                        number of conditions
%                        {3}: [trial x stimulus], condition of each trial

blanks = variables{2};
conditions = variables{3};
notBlanks = true(size(respTraining));
notBlanks(:,blanks) = false;
if nargin < 4
    testIndices = find(isnan(respTraining));
end
valid = find(~isnan(conditions) & notBlanks);
ind = intersect(testIndices, valid);

conds = unique(conditions(ind));
prediction = NaN(length(ind),1);
for a = 1:length(conds)
    k = conditions(ind)==conds(a);
    prediction(k) = nanmean(respTraining(conditions==conds(a) & notBlanks));
end
error = NaN(length(testIndices), 1);
error(ismember(testIndices, valid)) = prediction - respTest(ind);