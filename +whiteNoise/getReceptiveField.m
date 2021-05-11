function [receptiveFields, runKernels, runWin, ...
    explainedVariance, explainedVariance_runOnly, ...
    explainedVariance_stimOnly, predictions, predictions_runOnly, time] = ...
    getReceptiveField(traces, traceTimes, ...
    stimFrames, stimTimes, RFtimesInFrames, ...
    runSpeed, runTime, runKrnlLimits, lambdas, crossFolds)

%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveField, runKernels, runWin, ...
%    explainedVariance, explainedVariance_runOnly, ...
%    explainedVariance_stimOnly, predictions, predictions_runOnly, time] = ...
%    GETRECEPTIVEFIELD(traces, ...
%    traceTimes, stimFrames, stimTimes, RFtimesInFrames, ...
%    runSpeed, runTime, runKrnlLimits, lambdas, crossFolds) calculates the linear RF of the neuron.
%
%   receptiveFields     [rows x cols x RFframes x RFtype x neuron]
%                       containing linear regression solution for x in Ax=B
%                       where A is stimulus [rows x cols x time] and B is
%                       calcium response, for each neuron and stimulus 
%                       model; ridge regression is performed on all data
%                       using the optimal lambda value found with
%                       cross-validation
%   runKernels          [t x neuron]; containing linear
%                       regression kernel fitting calcium based on running
%                       speed
%   runWin              [1 x t]; time of run kernel relative to neural
%                       response
%   explainedVariance   [neuron x crossFold x lambdaStim x lambdaRun], each entry:
%                       explained variance for fitted RF and run kernel for
%                       each neuron, cross val. fold, and combination of
%                       lambdas
%   explainedVariance_runOnly   [neuron x crossFold x lambdaStim x lambdaRun], 
%                       each entry: explained variance by run kernel only, 
%                       for each neuron, cross val. fold, and combination of lambdas
%   explainedVariance_stimOnly   [neuron x crossFold x lambdaStim x lambdaRun], 
%                       each entry: explained variance by stimulus kernel only, 
%                       for each neuron, cross val. fold, and combination of lambdas
%   predictions         [t x neuron], each column contains
%                       prediction based on RF and run kernel for test 
%                       responses of specific neuron (using optimal
%                       lambdas)
%   predictions_runOnly [t x neuron], each column contains
%                       prediction based on run kernel only for test 
%                       responses of specific neuron (using optimal
%                       lambdas)
%   time                [t x 1]; time points for predictions
%
%   traces              [trTime x neuron]; calcium traces of neurons
%   traceTimes          [trTime x 1]; sample times of calcium traces
%   stimFrames          [time x rows x cols]; noise stimulus
%   stimTimes           [time x 1]; times of stimulus frames
%   RFtimesInFrames     [1 x RFframes]; frames of receptive field relative
%                       to stimulus frames
%   runSpeed            [rTime x 1]; running speed
%   runTime             [rTime x 1]; time points of running speed
%   runKrnlLimits       [1 x 2]; time limits of running kernel in s
%   lambdas             [1 x lambda]; values of lambda
%   crossFolds          ind; number of cross val. folds

if ~iscell(lambdas)
    lambdas = {lambdas, lambdas};
end

% generate toeplitz matrix for stimuli: [time x pixels]
% each row holds all pixels at current and previous time points:
% [[all pixels at t=0], [all pixels at t=-1], ...]
% each column is time series of that particular pixel

% find time gaps in stimulus presentation (usually when same visual noise
% stimulus was repeated several times)
stimBin = median(diff(stimTimes));
indGap = find(diff(stimTimes) > 2 * stimBin);
time = stimTimes;
% fill gaps with zeros in stimulus matrix
for g = 1:length(indGap)
    add = round(diff(stimTimes(indGap(g) + [0 1])) ./ stimBin);
    stimFrames = [stimFrames(1:indGap(g),:,:); ...
        zeros(add, size(stimFrames,2), size(stimFrames,3)); ...
        stimFrames(indGap(g)+1:end,:,:)];
    time = [time(1:indGap(g)); ...
        time(indGap(g)) + (1:add)' .* stimBin; ...
        time(indGap(g)+1:end)];
end
% reshape stimulus frames to [time x px]; this represents a single
% "stimulus block", i.e. the pixels to estimate a single time point of the
% receptive field
stim = reshape(stimFrames, size(stimFrames,1), []);
% now concatinate time shifted stimulus blocks; for each time point there
% is a stimulus block for lag=0, another for lag=-1, another for lag=-2,...
st = [];
for t = 1:length(RFtimesInFrames)
    st = [st, ...
        [zeros(max(0,RFtimesInFrames(1)-1+t), size(stim,2)); ...
        stim(max(1,2-RFtimesInFrames(1)-t) : end-RFtimesInFrames(1)-t+1, :)]];
end
stim = st;
clear st

% generate Toeplitz matrix for running speed
% filter running speed
runSpeed = medfilt1(runSpeed, 5);
% resample running speed at time points of neural traces
runBin = median(diff(runTime));
numBins = round(stimBin / runBin);
runSpeed = smooth(runSpeed, numBins);
runSpeed = interp1(runTime, runSpeed, time, 'pchip');
[runToepl, ~, runWin] =  krnl.getToeplitz(time, [], [], {runSpeed}, ...
    {runKrnlLimits}, true);
runWin = runWin{1};
% z-score
runToepl = (runToepl - mean(runToepl(:))) ./ std(runToepl);

% get neural response
traceBin = median(diff(traceTimes));
numBins = round(stimBin / traceBin);
traces = smoothdata(traces, 1, 'movmean', numBins, 'omitnan');
traces = interp1(traceTimes, traces, time);
% z-score neural response
zTraces = (traces - nanmean(traces,1)) ./ nanstd(traces,0,1);

% delete stim frames for which all neurons have NaN
ind = all(isnan(zTraces),2);
stim(ind,:) = [];
zTraces(ind,:) = [];
runToepl(ind,:) = [];
time(ind) = [];
% if NaN values < 5% in a neuron, exchange NaNs for 0
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.05;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
% skip neurons that have only NaN values
valid = ~all(isnan(zTraces),1)';

% duplicate stimulus matrix to predict ON part (1st half) and OFF
% part (2nd half)
s = stim;
s(stim < 0) = 0;
stim2 = s;
s = stim;
s(stim > 0) = 0;
stim2 = [stim2, s];
stim2 = (stim2 - nanmean(stim2(:))) ./ nanstd(stim2(:)); % normalise each column of stimulus matrix
clear s

% scale lamdas according to number of samples and number of predictors
lamStim = sqrt(lambdas{1} .* size(stim,1) .* size(stim,2));
lamRun = sqrt(lambdas{2} .* size(stim,1) .* size(runToepl,2));

% construct spatial smoothing lambda matrix
lamMatrix_stim = krnl.makeLambdaMatrix([size(stimFrames,2), size(stimFrames,3), ...
    length(RFtimesInFrames)], [1 1 0]);
lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);
lamMatrix_run = krnl.makeLambdaMatrix(size(runToepl,2), 1);

nPerFold = ceil(size(stim,1) / crossFolds);

runKernels = NaN(length(runWin), size(traces,2));
explainedVariance = NaN(size(traces,2), crossFolds, length(lamStim), length(lamRun));
predictions = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim), length(lamRun));
explainedVariance_runOnly = NaN(size(traces,2), crossFolds, length(lamStim), length(lamRun));
predictions_runOnly = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim), length(lamRun));
explainedVariance_stimOnly = NaN(size(traces,2), crossFolds, length(lamStim), length(lamRun));

residuals = zTraces;

% get variances explained by stimulus alone and by stimulus and running
fprintf('  Folds (of %d) to get expl. var. of RF: ', crossFolds)
for fold = 1:crossFolds
    fprintf('%d ',fold)
    ind = (1:nPerFold) + (fold-1)*nPerFold;
    ind(ind > size(zTraces,1)) = [];
    j = true(size(zTraces,1),1);
    j(ind) = false;
    
    y_train = gpuArray(padarray(residuals(j,valid), size(lamMatrix_stim,1), 'post'));
    y_test = residuals(~j,valid);
    x_train = stim2(j,:);
    x_test = stim2(~j,:);

    y2_train = padarray(y_train, size(lamMatrix_run,1), 'post');
    x2_train = [x_train, runToepl(j,:)];
    x2_test = [x_test, runToepl(~j,:)];

    for lamS = 1:length(lamStim)
        lms = lamMatrix_stim .* lamStim(lamS);
        
        for lamR = 1:length(lamRun)
            lmr = lamMatrix_run .* lamRun(lamR);
            A = gpuArray([x2_train; ...
                [[lms, zeros(size(lms,1), size(lmr,2))]; ...
                [zeros(size(lmr,1), size(lms,2)), lmr]]]);
            
            B = gather(A \ y2_train);
            pred = x2_test * B; % get prediction
            predictions(1:sum(~j), fold, valid, lamS, lamR) = pred;
            explainedVariance(valid, fold, lamS, lamR) = 1 - ...
                sum((y_test - pred) .^ 2,1) ./ ...
                sum((y_test - mean(zTraces(j, valid),1)) .^2, 1);
            
            pred = x_test * B(1:size(stim2,2),:);
            explainedVariance_stimOnly(valid, fold, lamS, lamR) = 1 - ...
                sum((y_test - pred) .^ 2,1) ./ ...
                sum((y_test - mean(zTraces(j, valid),1)) .^ 2,1);
            
            pred = runToepl(~j,:) * B(size(stim2,2)+1 : end,:);
            predictions_runOnly(1:sum(~j), fold, valid, lamS, lamR) = pred;
            explainedVariance_runOnly(valid, fold, lamS, lamR) = 1 - ...
                sum((y_test - pred) .^ 2,1) ./ ...
                sum((y_test - mean(zTraces(j, valid),1)) .^ 2,1);
        end
    end
end
fprintf('\n')

% determine RFs using all data and optimal lambdas
receptiveFields = NaN(size(stim2,2), size(traces,2));

[evRun, bestStimLams] = max(mean(explainedVariance, 2), [], 3);
evRun = permute(evRun, [1 4 2 3]); % [neuron x 1/lambdasRun]
bestStimLams = permute(bestStimLams, [1 4 2 3]); % [neuron x 1/lambdasRun]
[~, bestRunLams] = max(evRun, [], 2); % [neuron x 1]
ind = sub2ind(size(bestStimLams), (1:size(bestStimLams,1))', bestRunLams);
bestStimLams = bestStimLams(ind); % [neuron x 1]
fprintf('  Optimal lambdas (of %d) to get RFs: ', length(lamStim))
for lamS = 1:length(lamStim)
    fprintf('%d ', lamS)
    ind1 = bestStimLams == lamS & valid;
    A = [stim2; lamMatrix_stim .* lamStim(lamS)];
    tr = padarray(residuals, size(lamMatrix_stim,1), 'post');
    
    for lamR = 1:length(lamRun)
        ind2 = ind1 & bestRunLams == lamR;
        A2 = [[A, padarray(runToepl, size(lamMatrix_stim,1), 'post')]; ...
            padarray(lamMatrix_run .* lamRun(lamR), [0 size(A,2)], 'pre')];
        tr2 = padarray(tr(:,ind2), size(lamMatrix_run,1), 'post');
        B = gather(gpuArray(A2) \ gpuArray(tr2));
        receptiveFields(:,ind2) = B(1:size(stim2,2),:); % get RF kernel
        runKernels(:,ind2) = B(size(stim2,2)+1:end,:);
        predictions(:,:,ind2,1,1) = predictions(:,:,ind2,lamS,lamR);
        predictions_runOnly(:,:,ind2,1,1) = predictions_runOnly(:,:,ind2,lamS,lamR);
    end
end
fprintf('\n')

receptiveFields = reshape(receptiveFields, size(stimFrames,2), ...
    size(stimFrames,3), length(RFtimesInFrames), 2, size(traces,2));

predictions = reshape(predictions(:,:,:,1,1), [], size(traces,2));
predictions = predictions(1:length(time),:);
predictions_runOnly = reshape(predictions_runOnly(:,:,:,1,1), [], size(traces,2));
predictions_runOnly = predictions_runOnly(1:length(time),:);