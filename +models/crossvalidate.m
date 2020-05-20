function errors = crossvalidate(model, neuralResponses, variables, mode)
% model             function that fits model and returns prediction error for
%                   test set (e.g. models.testFitlmModel)
% neuralResponses   {neurons x 1}; each entry: [trials x stimuli]
% variables         {neurons x numVars}; each entry: [trials x stimuli];
%                   if size(variables,1)==1, same data are used for each
%                   neurons
% mode              String, either 'leaveOneOut' or 'stimSets'; if
%                   'leaveOneOut' model is trained on all but one trial, if
%                   'stimSets' test set consists of one random trial from
%                   each stimulus

% errors            {neurons x 1}; each entry: [trials x stimuli];
%                   difference between prediction (when this trial was in
%                   test set) and data

errors = cell(size(neuralResponses));
for iNeuron = 1:size(neuralResponses,1)
    resp = neuralResponses{iNeuron};
    [numTrials, numStim] = size(resp);
    switch mode
        case 'leaveOneOut'
            indsPerSet = num2cell((1:numel(resp))');
        case 'stimSets'
            trials = zeros(size(resp));
            for k = 1:numStim
                trials(:,k) = randperm(numTrials);
            end
            indsPerSet = cell(numTrials,1);
            for k = 1:numTrials
                indsPerSet{k} = find(trials == k);
            end
    end
    err = cell(size(indsPerSet));
%     parfor k = 1:length(indsPerSet)
    for k = 1:length(indsPerSet) % use 'for' if Parallel computing toolbox not installed
        respTraining = resp;
        respTraining(indsPerSet{k}) = NaN;
        respTest = NaN(size(resp));
        respTest(indsPerSet{k}) = resp(indsPerSet{k});
        err{k} = model(respTraining, respTest, variables, indsPerSet{k}); %#ok<PFBNS>
    end
    errs = NaN(size(resp));
    for k = 1:length(indsPerSet)
        errs(indsPerSet{k}) = err{k};
    end
    errors{iNeuron} = errs;
end