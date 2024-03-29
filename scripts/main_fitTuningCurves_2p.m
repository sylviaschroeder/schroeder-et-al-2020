%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish';
folderTools = 'C:\STORAGE\workspaces';
folderThisRepo = 'C:\dev\workspace\schroeder-et-al-2020';
folderResults = 'C:\dev\test';

%% Parameters
dataset = 'boutons';
% dataset = 'sc neurons 2p';

% smoothing (low-pass filter) of pupil size before
smoothStd = 0.25; %in sec
% parameters to separate trials into low and high arousal
threshPerc = 50;
minDurPerTrial = .5;
labels = {'small pupil','large pupil'};

% curve fitting
degrees = 1:360;

doSave = 1;
doPlot = 1;

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Fit tuning curves to low and high arousal responses
% fit tuning curves using different constraints: define which curve
% parameters have to have the same values under both conditions;
% parameters: (1) pref. dir., (2) ampl. at pref. dir, (3) direction index,
%             (4) offset of tuning curve, (5) tuning width

% only parameter sets (4) and (5) were used for further analyses; the other
% sets were only used to compare explained variances
parameterSets = {[], 1, [1 5], [1 3 5], 'constant'};
fileNames = {'tuning_nothingFixedAcrossConditions.mat', ...
    'tuning_prefDirFixed.mat', ...
    'tuning_prefDirSigmaFixed.mat', ...
    'tuning_prefDirSigmaDIFixed.mat', ...
    'tuning_constantFit.mat'};

subjects = dir(fullfile(folderBase, dataset, 'SS*'));

for iPars = length(parameterSets):-1:1 % results from parameterSets 
    fixedPars = parameterSets{iPars};
    if strcmp(parameterSets, 'constant')
        nPars = 1;
    else
        nPars = 5;
    end
    
    for subj = 1:length(subjects)
        name = subjects(subj).name;
        dates = dir(fullfile(folderBase, dataset, name, '2*'));
        for dt = 1:length(dates)
            date = dates(dt).name;
            fprintf('Dataset: %s %s\n', name, date);
            folder = fullfile(folderBase, dataset, name, date, '001');
            
            fResults = fullfile(folderResults, dataset, name, date);
            if ~isfolder(fResults)
                mkdir(fResults)
            end
            
            % load data
            time = readNPY(fullfile(folder, '_ss_2pCalcium.timestamps.npy'));
            kernels = readNPY(fullfile(folder, '_ss_2pRois._ss_gratingsKernels.npy'));
            isGad = readNPY(fullfile(folder, '_ss_2pRois.isGad.npy'));
            pupilSize = readNPY(fullfile(folder, 'eye.diameter.npy'));
            pupilTime = readNPY(fullfile(folder, 'eye.timestamps.npy'));
            stimIntervals = readNPY(fullfile(folder, '_ss_grating.intervals.npy'));
            stimSequence = readNPY(fullfile(folder, '_ss_grating._ss_gratingID.npy'));
            directions = readNPY(fullfile(folder, '_ss_gratingID.directions.npy'));
            amplitudes = readNPY(fullfile(folder, '_ss_gratingTrials.amplitudes.npy'));
            blanks = isnan(directions);
            directions(blanks) = [];
            timeBin = median(diff(time));
            
            isSuppressed = NaN(length(isGad),1);
            crossValExplVar = NaN(length(isGad),1);
            pars_smallPupil = NaN(length(isGad), nPars);
            pars_largePupil = NaN(length(isGad), nPars);
            curves_smallPupil = NaN(length(isGad), 360);
            curves_largePupil = NaN(length(isGad), 360);
            
            stimMatrix = exp.buildStimMatrix(stimSequence, stimIntervals, time);
            stimDurFrames = round(mean(diff(stimIntervals,1,2)) / median(diff(time)));
            
            % determine condition (small/large pupil) for each trial
            % threshold between small and large pupil
            threshold = prctile(pupilSize,threshPerc);
            ind = isnan(pupilSize);
            timeEdges = [time - 0.5*timeBin; time(end)+0.5*timeBin];
            indInterp = histcounts(pupilTime(ind), timeEdges) > 0;
            warning off
            pupilSize = interp1(pupilTime, pupilSize, time, 'pchip');
            warning on
            pupilSize = double(pupilSize > threshold);
            pupilSize(indInterp) = NaN;
            pupilSize = exp.getTracesPerStimulus(pupilSize, stimMatrix, [1 0]);
            pupilSize = squeeze(nanmean(pupilSize,4)); % [stimuli x repetition]
            pupilSize = double(pupilSize >= minDurPerTrial) + 1;
            nonVisualNoBlanks = pupilSize;
            nonVisualNoBlanks(blanks,:) = [];
            
            largePupil = pupilSize == 2;
            writeNPY(largePupil, fullfile(fResults, '_ss_gratingTrials.largePupil.npy'))
            
            fprintf('  Cell (of %d):', length(isGad));
            for iCell = 1:length(isGad)
                fprintf(' %d', iCell)
                resp = squeeze(amplitudes(iCell,:,:));
                if all(isnan(resp(:)))
                    continue
                end
                
                kernelSign = sign(sum(kernels(iCell,1:stimDurFrames)));
                medianResp = nanmedian(resp,2);
                [~,maxInd] = max(abs(medianResp));
                respSign = sign(medianResp(maxInd));
                isSuppressed(iCell) = -kernelSign * respSign;
                % make sure that kernel during stimulus presentation is always
                % positive
                resp = resp * kernelSign;
                % for fitting, all responses should be positive, so that the
                % preferred direction of suppressed-by-contrast cells will be
                % at the most negative response
                respFit = resp * respSign * kernelSign;
                
                respNoBlanks = respFit;
                respNoBlanks(blanks,:) = [];
                curves = NaN(length(degrees),2);
                if strcmp(fixedPars, 'constant')
                    % crossvalidate
                    errors = models.crossvalidate( @models.fitConstant, ...
                        {respNoBlanks'}, {directions, [], ...
                        nonVisualNoBlanks'}, 'stimSets');
                    
                    parameters = NaN(1, 2);
                    predictions = NaN(size(respFit));
                    for c = 1:2
                        ind = nonVisualNoBlanks == c;
                        parameters(c) = nanmean(respNoBlanks(ind));
                        predictions(ind) = parameters(c);
                        curves(:,c) = parameters(c);
                    end
                    if kernelSign * respSign < 0
                        parameters = -parameters;
                        predictions = -predictions;
                        curves = -curves;
                    end
                else
                    % crossvalidate
                    errors = models.crossvalidate( ...
                        @gratings.fitTuningCurveConditions_forCrossVal, ...
                        {respNoBlanks'}, {directions, false(length(directions),1), ...
                        nonVisualNoBlanks', fixedPars}, 'stimSets');
                    
                    [parameters, predictions] = gratings.fitTuningCurveConditions( ...
                        respFit, directions, blanks, pupilSize, fixedPars);
                    if kernelSign * respSign < 0
                        parameters(2,:) = -parameters(2,:);
                        parameters(4,:) = -parameters(4,:);
                        predictions = -predictions;
                    end
                    for c = 1:2
                        curves(:,c) = gratings.orituneWrappedConditions(parameters(:,c),degrees);
                    end
                end
                
                errors{1} = errors{1}';
                ind = ~isnan(errors{1}) & ~isnan(respNoBlanks);
                crossValExplVar(iCell) = 1 - sum(errors{1}(ind).^2) / ...
                    sum((respNoBlanks(ind)-mean(respNoBlanks(ind))).^2);
                pars_smallPupil(iCell,1:size(parameters,1)) = parameters(:,1);
                pars_largePupil(iCell,1:size(parameters,1)) = parameters(:,2);
                curves_smallPupil(iCell,:) = curves(:,1);
                curves_largePupil(iCell,:) = curves(:,2);
            end
            writeNPY(isSuppressed, fullfile(fResults, '_ss_tuning.isSuppressed.npy'))
            writeNPY(crossValExplVar, fullfile(fResults, '_ss_tuning.explVar.npy'))
            writeNPY(pars_smallPupil, fullfile(fResults, '_ss_tuning.parametersSmall.npy'))
            writeNPY(pars_largePupil, fullfile(fResults, '_ss_tuning.parametersLarge.npy'))
            writeNPY(curves_smallPupil, fullfile(fResults, '_ss_tuning.curvesSmall.npy'))
            writeNPY(curves_largePupil, fullfile(fResults, '_ss_tuning.curvesLarge.npy'))
            fprintf('\n')
        end
    end
end
