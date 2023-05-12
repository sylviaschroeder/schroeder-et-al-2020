%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish\arousal_NYP_matlab';
folderTools = 'C:\STORAGE\workspaces';
folderThisRepo = 'C:\dev\workspace\schroeder-et-al-2020';
folderResults = 'C:\dev\test';

%% Parameters
dataset = 'boutons';
% dataset = 'sc neurons 2p';

% for correcting baseline drifts of calcium traces at start of experiments
driftWin = 20; % in s, window to test whether baseline is higher than normal
driftThresh = 1.5; % in std, threshold for drift
correctWin = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
% kernels and RFs simultaneously
lambdasStim = logspace(-4, 1, 6);
lambdasRun = logspace(0, 6, 7);
RFlimits = [0.2 0.4];
crossFolds = 10;

% parameters for running speed as predictor
runKrnlLimits = [-5 5];

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Fit RFs and get cross-validated explained variance

subjects = dir(fullfile(folderBase, dataset, 'SS*'));
for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folderBase, dataset, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        folder = fullfile(folderBase, dataset, name, date, '001');
        
        fResults = fullfile(folderResults, dataset, name, date);
        if ~isfolder(fResults)
            mkdir(fResults)
        end
        
        % load data
        traces = readNPY(fullfile(folder, '_ss_2pCalcium.dff.npy'));
        time = readNPY(fullfile(folder, '_ss_2pCalcium.timestamps.npy'));
        delays = readNPY(fullfile(folder, '_ss_2pPlanes.delay.npy'));
        planes = readNPY(fullfile(folder, '_ss_2pRois._ss_2pPlanes.npy'));
        runSpeed = readNPY(fullfile(folder, '_ss_running.speed.npy'));
        runTime = readNPY(fullfile(folder, '_ss_running.timestamps.npy'));
        stimTimes = readNPY(fullfile(folder, '_ss_sparseNoise.times.npy'));
        stimMaps = readNPY(fullfile(folder, '_ss_sparseNoiseID.map.npy'));
        stimSeq = readNPY(fullfile(folder, '_ss_sparseNoise._ss_sparseNoiseID.npy'));
        stimPos = readNPY(fullfile(folder, '_ss_sparseNoiseArea.edges.npy'));
        
        % interpolate calcium traces to align all to same time
        timeBin = median(diff(time));
        for d = 2:length(delays)
            indUnits = find(planes == d);
            for n = indUnits'
                if all(isnan(traces(:,n)))
                    continue
                end
                nanInd1 = isnan(traces(:,n));
                traces(:,n) = interp1(time(~nanInd1) + delays(d), ...
                    traces(~nanInd1,n), time, 'pchip');
                nanInd2 = histcounts(time(nanInd1) + delays(d), time) > 0;
                traces(nanInd2,n) = NaN;
            end
        end
        
        % remove strong baseline decay at start of experiment in cells that
        % show it
        indUnits = find(nanmean(traces(1:round(driftWin / timeBin),:),1) > ...
            nanmean(traces,1) + driftThresh .* nanstd(traces,0,1));
        ind = round(correctWin / timeBin);
        for iUnit = 1:length(indUnits)
            y = traces(:, indUnits(iUnit));
            y = fillmissing(y, 'linear');
            % fit double exponential to start of trace
            f = fit((1:length(y))', y, ...
                @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
                'Lower', [0 0 0 0 0], ...
                'Upper', [max(y) max(y) 500 max(y) 500], ...
                'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
            % remove fit
            traces(:, indUnits(iUnit)) = y - f(1 : size(traces,1)) + f.a;
        end
        
        % more stimulus information
        stimFrames = stimMaps(stimSeq,:,:);
        stimFrameDur = median(diff(stimTimes));
        RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
            ceil(RFlimits(2) / stimFrameDur);
        clear stimMaps stimSeq
        
        % only consider left (contralateral) side of stimulus
        numCols = size(stimFrames,3);
        stimFrames = stimFrames(:,:,1:min(34, numCols));
        rfPos = stimPos;
        rfPos(2) = rfPos(1) + diff(stimPos([1 2])) / numCols * min(34, numCols);
        
        % map RF
        [rFields, runKernels, runWin, ev, ev_run, ev_stim] = ...
            whiteNoise.getReceptiveField( ...
            traces, time, stimFrames, stimTimes, ...
            RFtimesInFrames, runSpeed, runTime, runKrnlLimits, ...
            {lambdasStim, lambdasRun}, crossFolds);
        
        v = squeeze(mean(ev,2)); % [neuron x lamStim x lamRun], average across cross-folds
        [maxEV, maxStimLam] = max(v,[],2);
        maxEV = squeeze(maxEV); % [neuron x lamRun];
        maxStimLam = squeeze(maxStimLam); % [neuron x lamRun];
        [maxEV, maxRunLam] = max(maxEV, [], 2); % [neuron x 1]
        indLam = sub2ind(size(maxStimLam), (1:size(maxStimLam,1))', maxRunLam);
        maxStimLam = maxStimLam(indLam); % [neuron x 1]
        
        vRun = squeeze(mean(ev_run,2)); % [neuron x lamStim x lamRun], average across cross-folds
        vStim = squeeze(mean(ev_stim,2)); % [neuron x lamStim x lamRun]
        inds = sub2ind(size(vRun), (1:size(vRun,1))', maxStimLam, maxRunLam);
        maxEVRun = vRun(inds); % [neuron x 1]
        maxEVStim = vStim(inds); % [neuron x 1]
        
        % test signficance of each RF
        [ev, ev_shift] = ...
            whiteNoise.receptiveFieldShiftTest( ...
            traces, time, stimFrames, stimTimes, ...
            RFtimesInFrames, runSpeed, runTime, ...
            runKernels, runWin, rFields, maxStimLam, 500);
        pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
        pvals(isnan(ev)) = NaN;
        
        writeNPY(permute(rFields, [5 1 2 3 4]), fullfile(fResults, '_ss_rf.maps.npy'))
        writeNPY(maxEV, fullfile(fResults, '_ss_rf.explVars.npy'))
        writeNPY(maxEVRun, fullfile(fResults, '_ss_rf.explVarsRunning.npy'))
        writeNPY(maxEVStim, fullfile(fResults, '_ss_rf.explVarsStim.npy'))
        writeNPY(maxRunLam, fullfile(fResults, '_ss_rf.lambdasRunning.npy'))
        writeNPY(maxStimLam, fullfile(fResults, '_ss_rf.lambdasStim.npy'))
        writeNPY(rfPos, fullfile(fResults, '_ss_rfDescr.edges.npy'))
        writeNPY(RFtimesInFrames * stimFrameDur, fullfile(fResults, '_ss_rfDescr.timestamps.npy'))
        writeNPY(runKernels, fullfile(fResults, '_ss_rfRunningKernels.dff.npy'))
        writeNPY(runWin', fullfile(fResults, '_ss_rfRunningKernels.timestamps.npy'))
        
        writeNPY(pvals, fullfile(fResults, '_ss_rf.pValues.npy'))
    end
end