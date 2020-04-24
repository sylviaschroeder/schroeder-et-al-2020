%% Folders
folderBase = '';
folderExp = fullfile(folderBase, 'sylvia', 'Subjects', '%s', '%s', '%s', 'alf');

%% Parameters
% smoothing (low-pass filter) of pupil size before
smoothStd = 0.25; %in sec
% parameters to separate trials into low and high arousal
threshPerc = 50;
minDurPerTrial = .5;
labels = {'small pupil','large pupil'};

doSave = 1;
doPlot = 1;

%% Fit tuning curves to low and high arousal responses
% fit tuning curves using different constraints: define which curve
% parameters have to have the same values under both conditions;
% parameters: (1) pref. dir., (2) ampl. at pref. dir, (3) direction index,
%             (4) offset of tuning curve, (5) tuning width
parameterSets = {[], [1], [1 5], [1 3 5], 'constant'};
fileNames = {'tuning_nothingFixedAcrossConditions.mat', ...
    'tuning_prefDirFixed.mat', ...
    'tuning_prefDirSigmaFixed.mat', ...
    'tuning_prefDirSigmaDIFixed.mat', ...
    'tuning_constantFit.mat'};

for iPars = [4 5 1:3] %length(parameterSets):-1:1
    fixedPars = parameterSets{iPars}; %[1 5]; %[1];
    x = 0:359;
    tuning = struct([]);
    % tuning.plane.isGad
    %             .isSuppressed
    %             .R2
    %             .cond.name
    %                  .cell.parameters
    %                  .cell.curve
    %                  .cell.directions
    %                  .cell.medians
    %                  .cell.meanAbsDev
    %                  .cell.prefDir
    %                  .cell.width
    %                  .cell.amplitude
    %                  .cell.maxMin
    %                  .cell.respAtOppDir
    %                  .cell.respAtOrthDir
    %                  .cell.blankMedian
    %                  .cell.blankMeanAbsDev
    
    for iExp = 1:length(results)
        fprintf('Dataset %d: %s %s exp.: %d\n', iExp, results(iExp).subject, ...
            results(iExp).date, results(iExp).expGratings);
        folder = [folderROIData filesep results(iExp).subject filesep ...
            results(iExp).date filesep num2str(results(iExp).expGratings)];
        fileStart = [results(iExp).date '_' num2str(results(iExp).expGratings) '_' ...
            results(iExp).subject];
        file = [fileStart '_2P_plane%03d_ROI.mat'];
        tuning(iExp).subject = results(iExp).subject;
        tuning(iExp).date = results(iExp).date;
        tuning(iExp).exp = results(iExp).expGratings;
        tuning(iExp).planes = results(iExp).planes;
        for iPlane = 1:length(results(iExp).plane)
            display(['  Plane ' num2str(results(iExp).planes(iPlane))])
            % load meta
            data=load(fullfile(folder, sprintf(file,results(iExp).planes(iPlane))));
            meta = data.meta;
            meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
                'cortexlab.net');
            tuning(iExp).plane(iPlane).cellIDs = results(iExp).plane(iPlane).cellIDs;
            if isfield(results(iExp).plane, 'isGad')
                tuning(iExp).plane(iPlane).isGad = results(iExp).plane(iPlane).isGad;
            end
            tuning(iExp).plane(iPlane).R2 = struct([]);
            tuning(iExp).plane(iPlane).isSuppressed = ...
                NaN(length(results(iExp).plane(iPlane).cellIDs),1);
            tuning(iExp).plane(iPlane).crossValExplVar = ...
                NaN(length(results(iExp).plane(iPlane).cellIDs),1);
            if iPlane == 1
                nonVisData = [];
                % load ball or pupil data
                switch nonVisualSignal
                    case 'running'
                        ballData = nonVis.getRunningSpeed(meta);
                        if ~isempty(ballData)
                            stdSamples = round(smoothStd / median(diff(ballData.t)));
                            convWindow = normpdf(-3*stdSamples:3*stdSamples, 0, stdSamples);
                            nonVisData = conv(ballData.total, convWindow, 'same');
                            nonVisData = nonVisData / median(diff(ballData.t)) / 53;
                            nonVisTime = ballData.t;
                        end
                        ylab = 'Running speed (cm/s)';
                    case 'pupil'
                        [pupilData, nonVisTime] = nonVis.loadPupilData(meta);
                        if ~isempty(pupilData)
                            nonVisTime(length(pupilData.x)+1:end) = [];
                            nonVisData = nonVis.getPupilDiam(pupilData);
                        end
                        threshold = prctile(nonVisData,threshPerc);
                        ylab = 'Pupil diameter (a.u.)';
                end
                if doPlot == 1
                    figure('Position',[1925 685 1915 420])
                    plot(nonVisTime, nonVisData)
                    hold on
                    plot(nonVisTime([1 end]), [1 1] * threshold)
                    xlabel('Time (s)')
                    ylabel(ylab)
                    title(sprintf('Dataset %d: %s %s exp.: %d\n', iExp, ...
                        results(iExp).subject, results(iExp).date, ...
                        results(iExp).expGratings),'Interpreter','none')
                    xlim(nonVisTime([1 end]))
                end
            end
            [~, stimSeq, stimMatrix, frameTimes] = ...
                ssLocal.getStimulusResponseInfo(meta);
            [directions, blanks] = gratings.getOrientations(stimSeq);
            stimEnd = round(results(iExp).plane(iPlane).stimDuration / ...
                median(diff(results(iExp).plane(iPlane).kernelTime)));
            ind = isnan(nonVisData);
            indInterp = hist(nonVisTime(ind), frameTimes) > 0;
            warning off
            nonVisual = interp1(nonVisTime, nonVisData, frameTimes, 'pchip')';
            nonVisual = double(nonVisual > threshold);
            nonVisual(indInterp) = NaN;
            nonVisual = ssLocal.getTracesPerStimulus(nonVisual, stimMatrix, [1 0]);
            nonVisual = squeeze(nanmean(nonVisual,4)); % [stimuli x repetition]
            nonVisual = double(nonVisual >= minDurPerTrial) + 1;
            nonVisualNoBlanks = nonVisual;
            nonVisualNoBlanks(blanks,:) = [];
            % conditions: 1: not running/small pupil, 2: running/large pupil
            tun = struct([]);
            tun(1).cond(1).name = labels{1};
            tun(1).cond(2).name = labels{2};
            fprintf('    Cell (of %d):', length(results(iExp).plane(iPlane).kernelFit));
            for iCell = 1:length(results(iExp).plane(iPlane).kernelFit)
                fprintf(' %d', iCell)
                resp = results(iExp).plane(iPlane).kernelFit(iCell).alphaEachTrial';
                if isempty(resp)
                    tun.cond(1).cell(iCell).curve = [];
                    tun.cond(2).cell(iCell).curve = [];
                    continue
                end
                
                kernelSign = sign(sum(results(iExp).plane(iPlane)...
                    .kernelFit(iCell).kernel(1:stimEnd)));
                medianResp = nanmedian(resp,2);
                [~,maxInd] = max(abs(medianResp));
                respSign = sign(medianResp(maxInd));
                tuning(iExp).plane(iPlane).isSuppressed(iCell) = -kernelSign * respSign;
                % make sure that kernel during stimulus presentation is always
                % positive
                resp = resp * kernelSign;
                % for fitting, all responses should be positive, so that the
                % preferred direction of suppressed-by-contrast cells will be
                % at the most negative response
                respFit = resp * respSign * kernelSign;
                
                respNoBlanks = respFit;
                respNoBlanks(blanks,:) = [];
                curves = NaN(length(x),2);
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
                    R2 = 1 - nansum((respFit(:)-predictions(:)).^2) / ...
                        nansum((respFit(:)-nanmean(respFit(:))).^2);
                else
                    % crossvalidate
                    errors = models.crossvalidate( ...
                        @gratings.fitTuningCurveConditions_forCrossVal, ...
                        {respNoBlanks'}, {directions, [], nonVisualNoBlanks', ...
                        fixedPars}, 'stimSets');
                    
                    [parameters, ~, predictions,R2] = gratings.fitTuningCurveConditions( ...
                        respFit, directions, blanks, nonVisual, fixedPars);
                    if kernelSign * respSign < 0
                        parameters(2,:) = -parameters(2,:);
                        parameters(4,:) = -parameters(4,:);
                        predictions = -predictions;
                    end
                    for c = 1:2
                        curves(:,c) = gratings.orituneWrappedConditions(parameters(:,c),x);
                    end
                end
                
                errors{1} = errors{1}';
                ind = ~isnan(errors{1}) & ~isnan(respNoBlanks);
                explVar = 1 - sum(errors{1}(ind).^2) / ...
                    sum((respNoBlanks(ind)-mean(respNoBlanks(ind))).^2);
                tuning(iExp).plane(iPlane).crossValExplVar(iCell) = explVar;
                
                baselinePerTrial = ssLocal.getTracesPerStimulus( ...
                    results(iExp).plane(iPlane).kernelFit(iCell).baseline, ...
                    stimMatrix, [0 0]);
                baselinePerTrial = squeeze(mean(baselinePerTrial,4));
                inds = 1:size(resp,1);
                inds(blanks) = [];
                R2ToBlank = 1 - nansum(reshape((resp(inds,:)-predictions(inds,:)).^2,[],1)) / ...
                    nansum(reshape((resp(inds,:)-nanmean(baselinePerTrial(blanks,:))).^2,[],1));
                tuning(iExp).plane(iPlane).R2(iCell).comparedToMean = R2;
                tuning(iExp).plane(iPlane).R2(iCell).comparedToBlankResp = R2ToBlank;
                resp = resp(directions(:,2),:);
                nv = nonVisual(directions(:,2),:);
                for c = 1:2
                    tun.cond(c).cell(iCell).parameters = parameters(:,c);
                    tun.cond(c).cell(iCell).curve = curves(:,c);
                    tun.cond(c).cell(iCell).directions = directions(:,1);
                    r = NaN(size(resp));
                    r(nv==c) = resp(nv==c);
                    tun.cond(c).cell(iCell).responses = r;
                    
                    tun.cond(c).cell(iCell).blankResponses = ...
                        baselinePerTrial(blanks,nonVisual(blanks,:)==c);
                end
            end
            fprintf('\n')
            for c = 1:2
                tuning(iExp).plane(iPlane).cond(c) = tun.cond(c);
            end
            if doSave == 1
                if ~exist(fullfile(folderResults, nonVisualSignal), 'dir')
                    mkdir(fullfile(folderResults, nonVisualSignal));
                end
                save(fullfile(folderResults, nonVisualSignal, ...
                    fileNames{iPars}), 'tuning', 'x')
            end
        end
    end
end

%% Compare tuning curve fits (which parameters fixed)
fileNames = {'tuning_nothingFixed.mat', ...
    'tuning_prefDirFixed.mat', ...
    'tuning_prefDirSigmaFixed.mat', ...
    'tuning_prefDirSigmaDIFixed.mat', ...
    'tuning_constantFit.mat'};
groups = {'-','D_p','D_p, sigma','D_p, sigma, DI','const'};

explVar = cell(size(fileNames));
for iPars = 1:length(fileNames)
    data = load(fullfile(folderResults, nonVisualSignal, fileNames{iPars}));
    tuning = data.tuning;
    for iExp = 1:length(tuning)
        for iPlane = 1:length(tuning(iExp).plane)
            explVar{iPars} = [explVar{iPars}; ...
                tuning(iExp).plane(iPlane).crossValExplVar];
        end
    end
end

explVar = cell2mat(explVar);
g = size(explVar,2);
bins = -.175:.05:1;
figure
for i = 1:g
    for j = 1:g
        subplot(g,g,(i-1)*g+j)
        if i == j
            hist(explVar(:,i), bins)
            xlim([-.2 1])
        else
            plot(explVar(:,j),explVar(:,i),'k.')
            axis([-.2 1 -.2 1])
            axis square
        end
        if i == 1
            title(groups{j})
        end
        if j == 1
            ylabel(groups{i})
        end
        if i == g && j == ceil(g/2)
            xlabel('Explained variance')
        end
    end
end
figure
hold on
for i = 1:g
    n = hist(explVar(:,i),bins);
    plot(bins,n,'LineWidth',2)
end
legend(groups{1:end})
xlabel('Explained variance')
ylabel('# Neurons')

anova1(explVar(:,1:end-1));

%% Add to tuning structure which neurons are not tuned
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed.mat'));
tuning = data.tuning;
x = data.x;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_constantFit.mat'));
tunC = data.tuning;

for iExp = 1:length(tuning)
    for iPlane = 1:length(tuning(iExp).plane)
        tuned = NaN(length(tuning(iExp).plane(iPlane).cellIDs), 1);
        tuned(~isnan(tuning(iExp).plane(iPlane).crossValExplVar)) = 0;
        tuned(tuning(iExp).plane(iPlane).crossValExplVar > ...
            tunC(iExp).plane(iPlane).crossValExplVar) = 1;
        tuning(iExp).plane(iPlane).isTuned = tuned;
        ind = find(tuned == 0);
        for k = 1:length(ind)
            for c = 1:2
                par = tunC(iExp).plane(iPlane).cond(c).cell(ind(k)).parameters;
                tuning(iExp).plane(iPlane).cond(c).cell(ind(k)).parameters = ...
                    par;
                tuning(iExp).plane(iPlane).cond(c).cell(ind(k)).curve = ...
                    ones(length(x),1) .* par;
            end
            tuning(iExp).plane(iPlane).crossValExplVar(ind(k)) = ...
                tunC(iExp).plane(iPlane).crossValExplVar(ind(k));
            tuning(iExp).plane(iPlane).R2(ind(k)) = ...
                tunC(iExp).plane(iPlane).R2(ind(k));
        end
    end
end
save(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'), 'tuning', 'x');

%% Add line fits
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;
x = data.x;

for iExp = 1:length(tuning)
    for iPlane = 1:length(tuning(iExp).plane)
        if isempty(tuning(iExp).plane(iPlane).cellIDs)
            continue
        end
        for iCell = 1:length(tuning(iExp).plane(iPlane).cond(1).cell)
            dat1 = tuning(iExp).plane(iPlane).cond(1).cell(iCell);
            if isempty(dat1.parameters)
                continue
            end
            dat2 = tuning(iExp).plane(iPlane).cond(2).cell(iCell);
            if tuning(iExp).plane(iPlane).isTuned(iCell) == 1
                lineFit.intercept = dat2.parameters(4) - dat1.parameters(4) * ...
                    dat2.parameters(2) / dat1.parameters(2);
                lineFit.slope = dat2.parameters(2) / dat1.parameters(2);
            else
                lineFit.intercept = dat2.parameters - dat1.parameters;
                lineFit.slope = NaN;
            end
            tuning(iExp).plane(iPlane).lineFit(iCell) = lineFit;
        end
    end
end
save(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'), 'tuning', 'x');

%% Create tuning null distribution (shuffle conditions)
draws = 200;
fixedPars = [1 3 5];
runFitTwice = 1;
data = load(fullfile(folderResults, nonVisualSignal, 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    null(iExp).subject = tuning(iExp).subject;
    null(iExp).date = tuning(iExp).date;
    null(iExp).exp = tuning(iExp).exp;
    null(iExp).planes = tuning(iExp).planes;
    for iPlane = 1:length(tuning(iExp).plane)
        if isempty(tuning(iExp).plane(iPlane).cond)
            continue
        end
        null(iExp).plane(iPlane).cellIDs = tuning(iExp).plane(iPlane).cellIDs;
        null(iExp).plane(iPlane).isSuppressed = tuning(iExp).plane(iPlane).isSuppressed;
        null(iExp).plane(iPlane).isTuned = tuning(iExp).plane(iPlane).isTuned;
        
        neuron = find(~isnan(tuning(iExp).plane(iPlane).isTuned), 1);
        conditions = NaN(size(tuning(iExp).plane(iPlane).cond(1).cell(neuron).responses));
        for c = 1:2
            ind = ~isnan(tuning(iExp).plane(iPlane).cond(c).cell(neuron).responses);
            conditions(ind) = c;
        end
        directions = tuning(iExp).plane(iPlane).cond(1).cell(neuron).directions;
        directions = [directions, (1:length(directions))'];
        
        fprintf('  Cell (of %d):', length(tuning(iExp).plane(iPlane).isTuned));
        permutations = NaN(numel(conditions), draws);
        for k = 1:draws
            permutations(:,k) = randperm(numel(conditions));
        end
        for iCell = 1:length(tuning(iExp).plane(iPlane).isTuned)
            fprintf(' %d', iCell)
            if isnan(tuning(iExp).plane(iPlane).isTuned(iCell))
                continue
            end
            resp = NaN(size(conditions));
            for c = 1:2
                resp = nansum(cat(3, resp, ...
                    tuning(iExp).plane(iPlane).cond(c).cell(iCell).responses), 3);
            end
            if tuning(iExp).plane(iPlane).isTuned(iCell) == 0
                parameters = NaN(draws, 2);
                for k = 1:draws
                    conds = reshape(conditions(permutations(:,k)), size(conditions));
                    for c = 1:2
                        ind = conds == c;
                        parameters(k,c) = nanmean(resp(ind));
                    end
                end
            else
                if tuning(iExp).plane(iPlane).isSuppressed(iCell) > 0
                    resp = -resp;
                end
                parameters = NaN(draws, 2, 5);
                for k = 1:draws
                    conds = reshape(conditions(permutations(:,k)), size(conditions));
                    pars = gratings.fitTuningCurveConditions(resp, ...
                        directions, [], conds, fixedPars, runFitTwice);
                    parameters(k,:,:) = permute(pars, [3 2 1]);
                end
                if tuning(iExp).plane(iPlane).isSuppressed(iCell) > 0
                    parameters(:,:,2) = -parameters(:,:,2);
                    parameters(:,:,4) = -parameters(:,:,4);
                end
            end
            
            for c = 1:2
                null(iExp).plane(iPlane).cond(c).cell(iCell).parameters = ...
                    squeeze(parameters(:,c,:));
            end
        end
    end
    fprintf('\n')
    save(fullfile(folderResults, nonVisualSignal, ...
        'nullTuning_prefDirSigmaDIFixed.mat'), 'null')
end

%% Plot distributions of modulation measures + null distributions
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'nullTuning_prefDirSigmaDIFixed.mat'));
null = data.null;

minis = [];
maxis = [];
stimMeans =[];
nullMinisFit = [];
nullMaxisFit = [];
nullStimMeans = [];
isSuppr = [];
isGad = [];
[draws,~] = cellfun(@size,{null(1).plane(1).cond(1).cell.parameters});
draws = max(draws);
for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    for iPlane = 1:length(tuning(iExp).plane)
        if isempty(tuning(iExp).plane(iPlane).cellIDs)
            continue
        end
        data = tuning(iExp).plane(iPlane);
        neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        
        for j = 1:length(neurons)
            iCell = neurons(j);
            means = NaN(1,2);
            prefs = NaN(1,2);
            nonprefs = NaN(1,2);
            nullMeans = NaN(1,2,draws);
            nullPrefs = NaN(1,2,draws);
            nullNonprefs = NaN(1,2,draws);
            for c = 1:2
                pars = data.cond(c).cell(iCell).parameters;
                curve = data.cond(c).cell(iCell).curve;
                if length(pars) == 1 % not tuned
                    means(c) = pars;
                else
                    means(c) = mean(curve);
                    oris = mod(pars(1) + [0 90 180], 360);
                    sr = gratings.orituneWrappedConditions(pars, oris);
                    prefs(c) = sr(1);
                    if sr(1)-sr(2)>0
                        [nonprefs(c), ind] = min(sr(2:3));
                    else % suppressed neurons
                        [nonprefs(c), ind] = max(sr(2:3));
                    end
                    ind = ind+1;
                end
                
                pars = null(iExp).plane(iPlane).cond(c).cell(iCell).parameters;
                if size(pars,2) == 1 % not tuned
                    nullMeans(1,c,:) = pars;
                else
                    sr = NaN(size(pars,1), 3);
                    curves = NaN(size(pars,1), length(degrees));
                    for p = 1:size(pars,1)
                        oris = mod(pars(p,1) + [0 90 180], 360);
                        sr(p,:) = gratings.orituneWrappedConditions(pars(p,:), oris);
                        curves(p,:) = gratings.orituneWrappedConditions(pars(p,:), degrees);
                    end
                    nullMeans(1,c,:) = mean(curves,2);
                    nullPrefs(1,c,:) = sr(:,1);
                    nullNonprefs(1,c,:) = sr(:,ind);
                end
            end
            
            minis(end+1,:) = nonprefs;
            maxis(end+1,:) = prefs;
            stimMeans(end+1,:) = means;
            nullMinisFit(end+1,:,:) = nullNonprefs;
            nullMaxisFit(end+1,:,:) = nullPrefs;
            nullStimMeans(end+1,:,:) = nullMeans;
            
            isSuppr(end+1,:) = data.isSuppressed(iCell);
            if isfield(data, 'isGad')
                isGad(end+1,:) = data.isGad(iCell);
            end
        end
    end
end

mods = abs(maxis - minis);
nullMods = abs(nullMaxisFit - nullMinisFit);

% invert sign of responses of suppressed cells
minis(isSuppr==1,:) = -minis(isSuppr==1,:);
maxis(isSuppr==1,:) = -maxis(isSuppr==1,:);
stimMeans(isSuppr==1,:) = -stimMeans(isSuppr==1,:);
nullMinisFit(isSuppr==1,:,:) = -nullMinisFit(isSuppr==1,:,:);
nullMaxisFit(isSuppr==1,:,:) = -nullMaxisFit(isSuppr==1,:,:);
nullStimMeans(isSuppr==1,:,:) = -nullStimMeans(isSuppr==1,:,:);

measures = {minis, maxis, mods, stimMeans};
nullMeasures = {nullMinisFit, nullMaxisFit, nullMods, nullStimMeans};
labels = {'Minimum', 'Maximum', 'Modulation depth', 'Mean'};
modFuns = @(a,b) (b-a)./(abs(a) + abs(b));
funLabels = 'Diff-index (large-small)/(large+small)';
binSizes = 0.05;

% Plot histograms and cumulative distributions of difference indices
for m = 1:length(measures)
    figure
    real = modFuns(measures{m}(:,1), measures{m}(:,2));
    pseudo = modFuns(squeeze(nullMeasures{m}(:,1,:)), ...
        squeeze(nullMeasures{m}(:,2,:)));
    confInt = prctile(pseudo, [2.5 97.5], 2);
    isSignf = real < confInt(:,1) | real > confInt(:,2);
    mini = round(min(real) / binSizes) * binSizes; %-.975;
    maxi = round(max(real) / binSizes) * binSizes; % .975;
    bins = mini:binSizes:maxi;
    n1 = hist(real(isSignf), bins);
    n2 = hist(real(~isSignf), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    maxi = max([-(mini-.5*binSizes) maxi+.5*binSizes]);
    mini = -maxi;
    xlim([mini maxi])
    r = nanmedian(abs(real));
    title(sprintf('%s (n=%d) (|b-a| ~ %d%% a)', labels{m}, ...
        sum(~isnan(real)), round(((1+r)/(1-r)-1)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    
    figure('Position', [1250 680 560 420])
    ind = ~isnan(real);
    x = sort(real(ind), 'ascend');
    y = (1:sum(ind)) ./ sum(ind);
    x = [mini; x; maxi];
    y = [0 y 1];
    ind = ~isnan(pseudo);
    xNull = sort(pseudo(ind), 'ascend');
    yNull = (1:sum(ind(:))) ./ sum(ind(:));
    xNull = [mini; xNull; maxi];
    yNull = [0 yNull 1];
    hold on
    plot(x, y, 'k', 'LineWidth', 2)
    plot(xNull, yNull, ':', 'Color', [.5 .5 .5], 'LineWidth', 2)
    xlim(x([1 end]))
    xlabel(funLabels)
    ylabel('Proportion of neurons')
    title(sprintf('%s (n=%d)', labels{m}, sum(~isnan(real))))
    legend('data','null distribution', 'Location', 'NorthWest')
    legend('boxoff')
end

% Plot scatter plots
% for s = 1:2
    m1 = 3; %(s-1)*2 + 1;
    m2 = 4; %(s-1)*2 + 2;
    real1 = modFuns(measures{m1}(:,1), measures{m1}(:,2));
    real2 = modFuns(measures{m2}(:,1), measures{m2}(:,2));
    mini = -.975; %round(min(real) / binSizes) * binSizes;
    maxi = .975; %round(max(real) / binSizes) * binSizes;
    figure
    scatter(real1, real2, 20, 'k', 'filled')
    alpha(.3)
    xlim([-1 1])
    ylim([-1 1])
    title(sprintf('n = %d', sum(~isnan(real1)&~isnan(real2))))
    xlabel(labels{m1})
    ylabel(labels{m2})
    axis square
% end

% Plot scatter plots for different cell types
real1(isnan(real1)) = 0;
real2(real2>.5) = .5;
real2(real2<-.5) = -.5;
inds = {isGad==1, isGad==-1, isSuppr==-1, isSuppr==1};
types = {'inhibitory','excitatory','enhanced','suppressed'};
cols = lines(5);
cols(3,:) = [];
for c = 1:length(inds)
    figure
    hold on
    plot([0 0],[-.5 .5],'k')
    plot([-1 1],[0 0],'k')
    scatter(real1(inds{c}), real2(inds{c}), 36, cols(c,:), 'filled', ...
        'MarkerFaceAlpha', .3)
    xlim([-1 1])
    ylim([-.5 .5])
    title(sprintf('%s (n = %d)', types{c}, sum(inds{c}&~isnan(real2))))
    xlabel(labels{m1})
    ylabel(labels{m2})
    axis square
end

%% Plot differences of tuning parameters across conditions
minEV = 0.25;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_nothingFixed.mat'));
tuning = data.tuning;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_constantFit.mat'));
tunC = data.tuning;
prefDirs = zeros(0,2);
DSIs = zeros(0,2);
OSIs = zeros(0,2);
DI = zeros(0,2); % DI as used to fit tuning curve (disregarding differences in offsets)
sigma = zeros(0,2); % as used to fit tuning curve
for iExp = 1:length(tuning)
    for iPlane = 1:length(tuning(iExp).plane)
        dat = tuning(iExp).plane(iPlane);
        for iCell = 1:length(dat.R2)
            if isnan(dat.crossValExplVar(iCell)) || ...
                    dat.crossValExplVar(iCell) < minEV || ...
                    dat.crossValExplVar(iCell) <= ...
                    tunC(iExp).plane(iPlane).crossValExplVar(iCell)
                continue
            end
            pars1 = dat.cond(1).cell(iCell).parameters;
            pars2 = dat.cond(2).cell(iCell).parameters;
            prefDirs(end+1,:) = [pars1(1), pars2(1)];
            DI(end+1,:) = [pars1(3), pars2(3)];
            sigma(end+1,:) = [pars1(5), pars2(5)];
            
            D1 = pars1(1) + [0 90 180];
            D2 = pars2(1) + [0 90 180];
            Rn1 = (1-pars1(3))/(1+pars1(3)) * pars1(2);
            Rn2 = (1-pars2(3))/(1+pars2(3)) * pars2(2);
            R1 = orituneWrapped([pars1(1:2); Rn1; pars1(4:5)], D1);
            R2 = orituneWrapped([pars2(1:2); Rn2; pars2(4:5)], D2);
            if dat.isSuppressed(iCell) == 1
                R1 = -R1;
                R2 = -R2;
            end
            mini = min([R1, R2]);
            if mini < 0
                R1 = R1 - mini;
                R2 = R2 - mini;
            end
            Rori1 = mean(R1([1 3]));
            Rori2 = mean(R2([1 3]));
            OSIs(end+1,:) = [(Rori1-R1(2))/(Rori1+R1(2)), ...
                (Rori2-R2(2))/(Rori2+R2(2))];
            DSIs(end+1,:) = [(R1(1)-R1(3))/(R1(1)+R1(3)), ...
                (R2(1)-R2(3))/(R2(1)+R2(3))];
        end
    end
end
figure
plot(prefDirs(:,1), prefDirs(:,2), 'k.')
axis([0 360 0 360])
hold on
plot([0 360],[0 360], 'k:')
plot([180 360],[0 180], 'k:')
plot([0 180],[180 360], 'k:')
axis square
set(gca, 'box', 'off', 'XTick', 0:90:360, 'YTick', 0:90:360)
xlabel(labels{1})
ylabel(labels{2})
title('Preferred direction')
figure
plot(DI(:,1), DI(:,2), 'k.')
axis([0 1 0 1])
hold on
plot([0 1],[0 1], 'k:')
axis square
set(gca, 'box', 'off')
xlabel(labels{1})
ylabel(labels{2})
title('DI (without offsets)')
figure
plot(sigma(:,1), sigma(:,2), 'k.')
axis([25 180 25 180])
hold on
plot([25 180],[25 180], 'k:')
axis square
set(gca, 'box', 'off')
xlabel(labels{1})
ylabel(labels{2})
title('Sigma')
figure
plot(OSIs(:,1), OSIs(:,2), 'k.')
axis([0 1 0 1])
hold on
plot([0 1],[0 1], 'k:')
axis square
set(gca, 'box', 'off', 'XTick', 0:.5:1, 'YTick', 0:.5:1)
xlabel(labels{1})
ylabel(labels{2})
title('OSI (with offsets)')
% figure
% plot(DSIs(:,1), DSIs(:,2), 'k.')
% axis([0 1 0 1])
% hold on
% plot([0 1],[0 1], 'k:')
% axis square
% set(gca, 'box', 'off')
% xlabel(labels{1})
% ylabel(labels{2})
% title('DSI (with offsets)')

% data = load(fullfile(folderResults, nonVisualSignal, ...
%     'tuning_prefDirSigmaFixed.mat'));
% tuning = data.tuning;
% DSIs = zeros(0,2);
% OSIs = zeros(0,2);
% for iExp = 1:length(tuning)
%     for iPlane = 1:length(tuning(iExp).plane)
%         dat = tuning(iExp).plane(iPlane);
%         for iCell = 1:length(dat.R2)
%             if isnan(dat.crossValExplVar(iCell)) || ...
%                     dat.crossValExplVar(iCell) < minEV || ...
%                     dat.crossValExplVar(iCell) <= ...
%                     tunC(iExp).plane(iPlane).crossValExplVar(iCell)
%                 continue
%             end
%             D = dat.cond(1).cell(iCell).parameters(1) + [0 90 180];
%             R1 = orituneWrapped(dat.cond(1).cell(iCell).parameters, D);
%             R2 = orituneWrapped(dat.cond(2).cell(iCell).parameters, D);
%             Rori1 = mean(R1([1 3]));
%             Rori2 = mean(R2([1 3]));
%             OSIs(end+1,:) = [(Rori1-R1(2))/(Rori1+R1(2)), ...
%                 (Rori2-R2(2))/(Rori2+R2(2))];
%             DSIs(end+1,:) = [(R1(1)-R1(3))/(R1(1)+R1(3)), ...
%                 (R2(1)-R2(3))/(R2(1)+R2(3))];
%         end
%     end
% end
% figure
% plot(OSIs(:,1), OSIs(:,2), 'k.')
% mini = 0;
% maxi = 1;
% axis([mini maxi mini maxi])
% hold on
% plot([mini maxi],[mini maxi], 'k:')
% axis square
% xlabel(labels{1})
% ylabel(labels{2})
% title('OSI')
% figure
% plot(DSIs(:,1), DSIs(:,2), 'k.')
% mini = 0;
% maxi = 1;
% axis([mini maxi mini maxi])
% hold on
% plot([mini maxi],[mini maxi], 'k:')
% axis square
% xlabel(labels{1})
% ylabel(labels{2})
% title('DSI')

%% Plot tuning curves and scatter plot with line fit for each neuron
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
% data = load(fullfile(folderResults, nonVisualSignal, ...
%     'tuning_prefDirSigmaDIFixed_bootstrap.mat'));
tuning = data.tuning;
x = data.x;
lin = {'-','-'};
cols = {'k', 'r'};
col2 = lines(1);
fPlots = fullfile(folderResults, ['plots_' nonVisualSignal], ...
    'tuningCurves');
if ~exist(fPlots, 'dir')
    mkdir(fPlots);
end
for iExp = 8:length(tuning)
    for iPlane = 1:length(tuning(iExp).plane)
        fPlots = fullfile(folderResults, ['plots_' nonVisualSignal], ...
            'tuningCurves', [tuning(iExp).subject '_' tuning(iExp).date ...
            '_' num2str(tuning(iExp).exp)], num2str(tuning(iExp).planes(iPlane)));
        if ~exist(fPlots, 'dir')
            mkdir(fPlots);
        end
        for iCell = 1:length(tuning(iExp).plane(iPlane).cellIDs)
            if isempty(tuning(iExp).plane(iPlane).cond(1).cell(iCell).parameters)
                continue
            end
            figure('Position', [680 680 1070 420])
%             figure('Position',[680 680 600 420])
            subplot(1,2,1)
            hold on
            h = [0 0];
            for c = 1:2
                h(c) = plot(x, tuning(iExp).plane(iPlane).cond(c) ...
                    .cell(iCell).curve, lin{c}, ...
                    'Color', cols{c}, 'LineWidth',2);
                errorbar(tuning(iExp).plane(iPlane).cond(c) ...
                    .cell(iCell).directions, ...
                    nanmean(tuning(iExp).plane(iPlane).cond(c) ...
                    .cell(iCell).responses,2), ...
                    nanstd(tuning(iExp).plane(iPlane).cond(c) ...
                    .cell(iCell).responses,0,2) ./ ...
                    sqrt(sum(~isnan(tuning(iExp).plane(iPlane).cond(c) ...
                    .cell(iCell).responses),2)), 'o', 'Color', cols{c})
                m = mean(tuning(iExp).plane(iPlane).cond(c) ...
                    .cell(iCell).blankResponses);
                s = std(tuning(iExp).plane(iPlane).cond(c) ...
                    .cell(iCell).blankResponses) / ...
                    length(tuning(iExp).plane(iPlane).cond(c) ...
                    .cell(iCell).blankResponses);
                fill(x([1 end end 1]), [[1 1].*(m+s), [1 1].*(m-s)], 'k', ...
                    'FaceColor', cols{c}, 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.3)
                plot(x([1 end]), [m m], '--', 'Color', cols{c}, 'LineWidth', 2)
            end
            set(gca,'XTick',0:90:360)
            legend(h, {tuning(iExp).plane(iPlane).cond.name})
            legend('boxoff')
            title(sprintf('EV = %.2f, R2_{stimResp} = %.2f, R2_{tuning} = %.2f', ...
                tuning(iExp).plane(iPlane).crossValExplVar(iCell), ...
                tuning(iExp).plane(iPlane).R2(iCell).comparedToBlankResp, ...
                tuning(iExp).plane(iPlane).R2(iCell).comparedToMean))
            xlim([-10 370])
            xlabel('Direction (in degrees)')
            ylabel('\DeltaF/F')
            
            subplot(1,2,2)
            resp1 = nanmean(tuning(iExp).plane(iPlane).cond(1) ...
                .cell(iCell).responses,2);
            resp2 = nanmean(tuning(iExp).plane(iPlane).cond(2) ...
                .cell(iCell).responses,2);
            if tuning(iExp).plane(iPlane).isSuppressed(iCell) == 1
                resp1 = -resp1;
                resp2 = -resp2;
            end
            maxiR = max(resp1);
            resp1 = resp1 ./ maxiR;
            resp2 = resp2 ./ maxiR;
            plot(resp1, resp2, 'ko')
            hold on
            intercept = tuning(iExp).plane(iPlane).lineFit(iCell).intercept / maxiR;
            slope = tuning(iExp).plane(iPlane).lineFit(iCell).slope;
            if isnan(slope)
                slope = 1;
            end
            if tuning(iExp).plane(iPlane).isSuppressed(iCell) == 1
                intercept = -intercept;
            end
            mini = min([0; resp1; resp2]);
            maxi = max([0; resp1; resp2]);
            rng = maxi-mini;
            mini = mini-.05*rng;
            maxi = maxi+.05*rng;
            plot([mini maxi],[mini maxi].*slope+intercept, ...
                'Color','k','LineWidth',2)
            plot([mini maxi],[mini maxi],'k:')
            axis([mini maxi mini maxi])
            axis square
            xlabel(tuning(iExp).plane(iPlane).cond(1).name)
            ylabel(tuning(iExp).plane(iPlane).cond(2).name)
            set(gca,'box','off')
            title(sprintf('y = %.2f x + %.2f', slope, intercept))
%             prcIntcpt = prctile(tuning(iExp).plane(iPlane).lineFit(iCell) ...
%                 .bsIntcpts, [2.5 97.5]) ./ maxiR;
%             if tuning(iExp).plane(iPlane).isSuppressed(iCell) == 1
%                 prcIntcpt = -prcIntcpt([2 1]);
%             end
%             prcSlope = prctile(tuning(iExp).plane(iPlane).lineFit(iCell) ...
%                 .bsSlps, [2.5 97.5]);
%             title(sprintf('y = %.2f [%.2f %.2f] x + %.2f [%.2f %.2f])', ...
%                 slope, prcSlope(1), prcSlope(2), intercept, prcIntcpt(1), ...
%                 prcIntcpt(2)))
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(fPlots, sprintf('tuningCurve%03d.jpg', iCell)), ...
                '-djpeg','-r0')
            close gcf
        end
    end
end

%% 

%% FROM HERE ON OLD

%% Bootstrap responses and fit tuning curves (old)
% Fix pref. dir., sigma, and DI across conditions (small/large pupil) as
% well as when doing bootstrap
draws = 250;

data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
x = data.x;

for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
            tuning(iExp).date, tuning(iExp).exp);
    for iPlane = 1:length(tuning(iExp).plane)
        display(['  Plane ' num2str(tuning(iExp).planes(iPlane))])
        fprintf('    Cell (of %d):', length(tuning(iExp).plane(iPlane).cond(1).cell));
        for iCell = 1:length(tuning(iExp).plane(iPlane).cond(1).cell)
            dat1 = tuning(iExp).plane(iPlane).cond(1).cell(iCell);
            if isempty(dat1.parameters)
                continue
            end
            fprintf(' %d', iCell)
            dat2 = tuning(iExp).plane(iPlane).cond(2).cell(iCell);
            resp = nansum(cat(3, dat1.responses, dat2.responses), 3);
            [numStim, numTrials] = size(resp);
            dirs = repmat(dat1.directions,1,numTrials);
            cond = ones(size(resp));
            cond(isnan(dat1.responses)) = 2;
            
%             sampling = randi(numTrials, numStim, numTrials, draws);
%             indStim = reshape(repmat((1:numStim)',1,numTrials), [], 1);
            sampling = randi(numel(resp), numel(resp), draws);
            bsPars = NaN(draws, 2, 2); % [draw x condition x parameter]
            if tuning(iExp).plane(iPlane).isTuned(iCell) == 1
                paramLimits = repmat(dat1.parameters',2,1);
                paramLimits(:,[2 4]) = repmat([-Inf;Inf],1,2);
                parfor d = 1:draws
                    %                 indTrial = reshape(sampling(:,:,d), [], 1);
                    %                 ind = sub2ind(size(sampling(:,:,d)), indStim, indTrial);
                    ind = sampling(:,d);
                    r = resp(ind);
                    c = cond(ind);
                    s = dirs(ind);
                    if length(unique(c)) < 2
                        continue
                    end
                    pars = fitoriConditions(s, r, c, [], paramLimits, 10);
                    pars = cumsum(reshape(pars, 5, []), 2)';
                    bsPars(d,:,:) = pars(:,[2 4]);
                end
            else
                parfor d = 1:draws
                    %                 indTrial = reshape(sampling(:,:,d), [], 1);
                    %                 ind = sub2ind(size(sampling(:,:,d)), indStim, indTrial);
                    ind = sampling(:,d);
                    r = resp(ind);
                    c = cond(ind);
                    if length(unique(c)) < 2
                        continue
                    end
                    bsPars(d,:,1) = [nanmean(r(c==1)), nanmean(r(c==2))];
                end
                bsPars(:,:,2) = [];
            end
            bootstrap(iExp).plane(iPlane).cond(1).cell(iCell).parameters = squeeze(bsPars(:,1,:));
            bootstrap(iExp).plane(iPlane).cond(2).cell(iCell).parameters = squeeze(bsPars(:,2,:));
        end
        fprintf('\n')
    end
end
save(fullfile(folderResults, nonVisualSignal, ...
    'bootstrap_prefDirSigmaDIFixed.mat'), 'bootstrap')

% for iExp = 1:3 %length(tuning)
%     for iPlane = 1:length(tuning(iExp).plane)
%         for iCell = 1:length(tuning(iExp).plane(iPlane).cond(1).cell)
%             dat1 = tuning(iExp).plane(iPlane).cond(1).cell(iCell);
%             if isempty(dat1.parameters)
%                 continue
%             end
%             dat2 = tuning(iExp).plane(iPlane).cond(2).cell(iCell);
%             bs1 = bootstrap(iExp).plane(iPlane).cond(1).cell(iCell).parameters;
%             bs2 = bootstrap(iExp).plane(iPlane).cond(2).cell(iCell).parameters;
%             if tuning(iExp).plane(iPlane).isTuned(iCell) == 1
%                 tuning(iExp).plane(iPlane).lineFit(iCell).ciIntcpts = ...
%                     prctile(bs2(:,2) - (bs1(:,2) .* ...
%                     bs2(:,1)) ./ bs1(:,1), [2.5 97.5]);
%                 tuning(iExp).plane(iPlane).lineFit(iCell).ciSlps = ...
%                     prctile(bs2(:,1) ./ bs1(:,1), [2.5 97.5]);
%             else
%                 tuning(iExp).plane(iPlane).lineFit(iCell).ciIntcpts = ...
%                     prctile(bs2-bs1, [2.5 97.5]);
%                 tuning(iExp).plane(iPlane).lineFit(iCell).ciSlps = NaN;
%             end
%         end
%     end
% end
% save(fullfile(folderResults, nonVisualSignal, ...
%     'tuning_prefDirSigmaDIFixed_bootstrap.mat'), 'tuning', 'x')

%% Plot other measures (reliability) across conditions
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;
FFs = zeros(0,2);
isTuned = zeros(0,1);
isSuppressed = zeros(0,1);
for iExp = 1:length(tuning)
    for iPlane = 1:length(tuning(iExp).plane)
        dat = tuning(iExp).plane(iPlane);
        for iCell = 1:length(dat.cellIDs)
            if isnan(dat.isTuned(iCell))
                continue
            end
            r1 = dat.cond(1).cell(iCell).responses;
            r2 = dat.cond(2).cell(iCell).responses;
            FFs(end+1,:) = [mean(nanvar(r1,[],2) ./ nanmean(r1,2)), ...
                mean(nanvar(r2,[],2) ./ nanmean(r2,2))];
            isTuned(end+1,:) = dat.isTuned(iCell);
            isSuppressed(end+1,:) = dat.isSuppressed(iCell);
        end
    end
end
FFs(isSuppressed==1,:) = -FFs(isSuppressed==1,:);
ind = any(FFs < 0, 2);
FFs(ind,:) = [];
isTuned(ind) = [];
isSuppressed(ind) = [];

figure
hold on
plot(FFs(isTuned==0,1), FFs(isTuned==0,2), '.', 'Color', [.8 .8 .8])
plot(FFs(isTuned==1,1), FFs(isTuned==1,2), 'k.')
xlabel(sprintf('Fano factor - %s', tuning(1).plane(1).cond(1).name))
ylabel(sprintf('Fano factor - %s', tuning(1).plane(1).cond(2).name))
mini = min(FFs(:));
maxi = max(FFs(:));
plot([mini maxi],[mini maxi],'r')
plot(median(FFs(:,1)), median(FFs(:,2)), 'ro', 'LineWidth', 2)
axis([0 10 0 10])
axis square
title(sprintf('Signed rank test: p = %.3f', signrank(FFs(:,1),FFs(:,2))))
legend('not tuned','tuned')

%% Plot bootstrapped tuning modulation measures (min, max, modulation)
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'bootstrap_prefDirSigmaDIFixed.mat'));
bootstrap = data.bootstrap;

minis = [];
maxis = [];
bootMinis = [];
bootMaxis = [];
isSuppr = [];
for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
            tuning(iExp).date, tuning(iExp).exp);
    for iPlane = 1:length(tuning(iExp).plane)
        neurons = find(tuning(iExp).plane(iPlane).isTuned == 1);
        for j = 1:length(neurons)
            iCell = neurons(j);
            mins = NaN(1,2);
            maxs = NaN(1,2);
            draws = size(bootstrap(iExp).plane(iPlane).cond(1) ...
                .cell(iCell).parameters,1);
            bMins = NaN(1,2,draws);
            bMaxs = NaN(1,2,draws);
            for c = 1:2
                pars = tuning(iExp).plane(iPlane).cond(c).cell(iCell).parameters;
                oris = mod(pars(1) + [0 90 180], 360);
                r = orituneWrappedConditions(pars, oris);
                maxs(c) = r(1);
                if r(1)-r(2)>=0
                    [mins(c), ind] = min(r(2:3));
                else % suppressed neurons
                    [mins(c), ind] = max(r(2:3));
                end
                ind = ind + 1;
                for k = 1:draws
                    p = pars;
                    p([2 4]) = bootstrap(iExp).plane(iPlane).cond(c) ...
                        .cell(iCell).parameters(k,:);
                    r = orituneWrappedConditions(p, oris);
                    bMins(1,c,k) = r(ind);
                    bMaxs(1,c,k) = r(1);
                end
            end
            minis(end+1,:) = mins;
            maxis(end+1,:) = maxs;
            bootMinis(end+1,:,:) = bMins;
            bootMaxis(end+1,:,:) = bMaxs;
            if maxs(1)>=mins(1)
                isSuppr(end+1,1) = 0;
            else
                isSuppr(end+1,1) = 1;
            end
        end
    end
end
modDepths = maxis - minis;
modDepths(isSuppr==1,:) = -modDepths(isSuppr==1,:);
bootModDepths = bootMaxis - bootMinis;
bootModDepths(isSuppr==1,:,:) = -bootModDepths(isSuppr==1,:,:);
diffsMinis = diff(minis,1,2);
diffsMinis(isSuppr==1) = -diffsMinis(isSuppr==1);
diffsBootMinis = squeeze(diff(bootMinis,1,2));
diffsBootMinis(isSuppr==1,:) = -diffsBootMinis(isSuppr==1,:);
diffsMaxis = diff(maxis,1,2);
diffsMaxis(isSuppr==1) = -diffsMaxis(isSuppr==1);
diffsBootMaxis = squeeze(diff(bootMaxis,1,2));
diffsBootMaxis(isSuppr==1,:) = -diffsBootMaxis(isSuppr==1,:);
diffsModDepths = diff(modDepths,1,2);
diffsBootModDepths = squeeze(diff(bootModDepths,1,2));

measures = {diffsMinis ./ sum(abs(minis),2), diffsMaxis ./ sum(abs(maxis),2), ...
    diffsModDepths ./ sum(modDepths,2)};
diffs = {diffsMinis, diffsMaxis, diffsModDepths};
boots = {diffsBootMinis, diffsBootMaxis, diffsBootModDepths};
labels = {'Minimum', 'Maximum', 'Modulation depth'};

binSize = 0.05;
mini = -0.975;
maxi = 0.975;
figure('Position', [40 600 1780 420])
for m = 1:length(measures)
    confInt = prctile(boots{m}, [2.5 97.5], 2);
    sign = confInt(:,1)>0 | confInt(:,2)<0;
    valid = diffs{m}>confInt(:,1) & diffs{m}<confInt(:,2);
    subplot(1,length(measures),m)
%     mini = round(min(measures{m}(valid)) / binSize) * binSize;
%     maxi = round(max(measures{m}(valid)) / binSize) * binSize;
    bins = mini:binSize:maxi;
    n1 = hist(measures{m}(valid & sign), bins);
    n2 = hist(measures{m}(valid & ~sign), bins);
    bar(bins, [n2',n1'], 'stacked')
    colormap([1 1 1;0 0 0])
    xlim([mini-.5*binSize maxi+.5*binSize])
    title(sprintf('%s (n=%d)', labels{m}, sum(valid)))
    if m == 1
        xlabel('Diff-index (large-small)/(large+small)')
        ylabel('#Neurons')
    end
    legend(sprintf('n=%d',sum(valid&~sign)), sprintf('n=%d',sum(valid&sign)))
end

%% Plot slopes and intercepts (population)
minExplVar = 0.2;

data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
% data = load(fullfile(folderResults, nonVisualSignal, ...
%     'tuning_prefDirSigmaDIFixed_bootstrap.mat'));
tuning = data.tuning;

% plot intercepts and slopes
% normalize all parameters so that max. response at small pupil/no running
% is 1
slopes = [];
intercepts = [];
% prctlSlopes = [];
% prctlIntercepts = [];
maxRespDiffs = [];
explVar = [];
isSuppr = [];
isGad = [];
isTuned = [];
for iExp = 1:length(tuning)
    for iPlane = 1:length(tuning(iExp).plane)
        data = tuning(iExp).plane(iPlane);
        neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        for j = 1:length(neurons)
            iCell = neurons(j);
            maxR = max(nanmean(data.cond(1).cell(iCell).responses, 2));
            slopes(end+1,1) = data.lineFit(iCell).slope;
            intercepts(end+1,1) = data.lineFit(iCell).intercept / maxR;
%             prctlSlopes(end+1,:) = data.lineFit(iCell).ciSlps;
%             prctlIntercepts(end+1,:) = data.lineFit(iCell).ciIntcpts;
            isSuppr(end+1,1) = data.isSuppressed(iCell);
            if isSuppr(end) == 1
                intercepts(end) = -intercepts(end);
%                 prctlIntercepts(end,:) = -prctlIntercepts(end,[2 1]);
            end
            isTuned(end+1,1) = data.isTuned(iCell);
            % diff. between responses at max. resp. in condition 1 (small
            % pupil/no running)
            if isTuned(end) == 1
                maxRespDiffs(end+1,1) = slopes(end) + intercepts(end) - 1;
            else
                maxRespDiffs(end+1,1) = intercepts(end);
            end
%             isGad(end+1,1) = data.isGad(iCell);
            explVar(end+1,1) = data.crossValExplVar(iCell);
        end
    end
end

% plot slopes against intercepts
ind = isTuned == 1 & explVar >= minExplVar;
xLimits = [-.62 .32];
yLimits = [0 2.9];
% xLimits = [-.2 .15];
% yLimits = [0.8 1.5];
figure
plot(intercepts(ind), slopes(ind) ,'k.')
hold on
plot([0 0], yLimits, 'k:')
plot(xLimits, [1 1], 'k:')
xlim(xLimits)
ylim(yLimits)
xlabel('Intercept')
ylabel('Slope')
title(sprintf('Tuned neurons (expl. var. >= %.2f)', minExplVar))
set(gca,'box','off')

% good = ind & (slopes>prctlSlopes(:,1) & slopes<prctlSlopes(:,2)) & ...
%     (intercepts>prctlIntercepts(:,1) & intercepts<prctlIntercepts(:,2));
%     
% onlyAdd = good & (prctlSlopes(:,1)<1 & prctlSlopes(:,2)>1) & ...
%     (prctlIntercepts(:,1)>0 | prctlIntercepts(:,2)<0);
% onlyMult = good & (prctlSlopes(:,1)>1 | prctlSlopes(:,2)<1) & ...
%     (prctlIntercepts(:,1)<0 & prctlIntercepts(:,2)>0);
% mixed = good & (prctlSlopes(:,1)>1 | prctlSlopes(:,2)<1) & ...
%     (prctlIntercepts(:,1)>0 | prctlIntercepts(:,2)<0);
% noChange = good & (prctlSlopes(:,1)<1 & prctlSlopes(:,2)>1) & ...
%     (prctlIntercepts(:,1)<0 & prctlIntercepts(:,2)>0);
% weirdos = (slopes<prctlSlopes(:,1) | slopes>prctlSlopes(:,2)) | ...
%     (intercepts<prctlIntercepts(:,1) | intercepts>prctlIntercepts(:,2));
% figure
% hold on
% plot(intercepts(noChange), slopes(noChange), '.', 'Color', [.5 .5 .5])
% plot(intercepts(mixed), slopes(mixed), 'g.')
% plot(intercepts(onlyMult), slopes(onlyMult), 'r.')
% plot(intercepts(onlyAdd), slopes(onlyAdd), 'b.')
% % plot(intercepts(weirdos), slopes(weirdos), 'c.')
% plot([0 0], yLimits, 'k:')
% plot(xLimits, [1 1], 'k:')
% xlim(xLimits)
% ylim(yLimits)
% xlabel('Intercept')
% ylabel('Slope')
% title(sprintf('Expl. var. >= %.2f', minExplVar))
% set(gca,'box','off')
% legend('no change','mixed','multiplicative','additive')

% plot slopes against intercepts for different cell classes
cols = lines(4);
ind = isTuned == 1 & explVar >= minExplVar;

gadPos = ind & isGad==1;
gadNeg = ind & isGad==-1;
figure
hold on
plot(intercepts(gadPos), slopes(gadPos), '.', 'Color', cols(1,:));
plot(intercepts(gadNeg), slopes(gadNeg), '.', 'Color', cols(2,:));
plot([0 0], yLimits, 'k:')
plot(xLimits, [1 1], 'k:')
xlim(xLimits)
ylim(yLimits)
xlabel('Intercept')
ylabel('Slope')
title(sprintf('Tuned neurons (expl. var. >= %.2f)', minExplVar))
set(gca,'box','off')
legend('Gad+','Gad-')

supprPos = ind & isSuppr==1;
supprNeg = ind & isSuppr==-1;
figure
hold on
plot(intercepts(supprNeg), slopes(supprNeg), '.', 'Color', cols(4,:));
plot(intercepts(supprPos), slopes(supprPos), '.', 'Color', cols(3,:));
plot([0 0], yLimits, 'k:')
plot(xLimits, [1 1], 'k:')
xlim(xLimits)
ylim(yLimits)
xlabel('Intercept')
ylabel('Slope')
title(sprintf('Tuned neurons (expl. var. >= %.2f)', minExplVar))
set(gca,'box','off')
legend('Suppressed-','Suppressed+')

% plot intercepts for untuned neurons
ind = isTuned == 0;
bins = -1.1:.05:1.6;
figure
hold on
hist(intercepts(ind), bins)
maxi = max(hist(intercepts(ind), bins));
plot([0 0], [0 maxi*1.2], 'k--')
plot(median(intercepts(ind)), maxi*1.1,'vk','MarkerFaceColor', 'k')
xlabel('Intercept')
ylabel('# Neurons')
title('Untuned neurons (line better than tuning curve)')
xlim(bins([1 end]))
ylim([0 maxi*1.2])
set(gca,'box','off')

ind = isTuned == 0;
gadPos = ind & isGad==1;
gadNeg = ind & isGad==-1;
figure
hold on
n1 = hist(intercepts(gadPos), bins);
n2 = hist(intercepts(gadNeg), bins);
plot(bins, n1, 'Color', cols(1,:), 'LineWidth', 2)
plot(bins, n2, 'Color', cols(2,:), 'LineWidth', 2)
maxi = max([n1,n2]);
plot([0 0], [0 maxi.*1.2], 'k--')
plot(median(intercepts(gadNeg)), maxi*1.1, 'v', 'Color', cols(2,:), ...
    'MarkerFaceColor', cols(2,:))
plot(median(intercepts(gadPos)), maxi*1.1, 'v', 'Color', cols(1,:), ...
    'MarkerFaceColor', cols(1,:))
ylim([0 1.2*maxi])
xlabel('Intercept')
ylabel('# Neurons')
title('Untuned neurons (line better than tuning curve)')
xlim(bins([1 end]))
set(gca,'box','off')
legend('Gad+','Gad-')

supprPos = ind & isSuppr==1;
supprNeg = ind & isSuppr==-1;
figure
hold on
n1 = hist(intercepts(supprNeg), bins);
n2 = hist(intercepts(supprPos), bins);
plot(bins, n1, 'Color', cols(4,:), 'LineWidth', 2)
plot(bins, n2, 'Color', cols(3,:), 'LineWidth', 2)
maxi = max([n1,n2]);
plot([0 0], [0 maxi.*1.2], 'k--')
plot(median(intercepts(supprNeg)), maxi*1.1, 'v', 'Color', cols(4,:), ...
    'MarkerFaceColor', cols(4,:))
plot(median(intercepts(supprPos)), maxi*1.1, 'v', 'Color', cols(3,:), ...
    'MarkerFaceColor', cols(3,:))
ylim([0 1.2*maxi])
xlabel('Intercept')
ylabel('# Neurons')
title('Untuned neurons (line better than tuning curve)')
xlim(bins([1 end]))
set(gca,'box','off')
legend('Suppressed-','Suppressed+')

% plot differences between maximal responses in both conditions
ind = (isTuned == 1 & explVar >= minExplVar) | isTuned == 0;
bins = -1.1:.05:1.1;
figure
hold on
hist(maxRespDiffs(ind), bins)
maxi = max(hist(maxRespDiffs(ind), bins));
plot([0 0], [0 maxi*1.2], 'k--')
plot(median(maxRespDiffs(ind)), maxi*1.1,'vk','MarkerFaceColor', 'k')
xlim(bins([1 end]))
ylim([0 maxi*1.2])
xlabel('Max. response diff.')
ylabel('# Neurons')
set(gca,'box','off')

gadPos = ind & isGad==1;
gadNeg = ind & isGad==-1;
n1 = hist(maxRespDiffs(gadPos), bins);
n2 = hist(maxRespDiffs(gadNeg), bins);
figure
hold on
plot(bins, n1, 'LineWidth', 2, 'Color', cols(1,:))
plot(bins, n2, 'LineWidth', 2, 'Color', cols(2,:))
maxi = max([n1,n2]);
plot([0 0], [0 maxi.*1.2], 'k--')
plot(median(maxRespDiffs(gadPos)), maxi*1.1, 'v', 'Color', cols(1,:), ...
    'MarkerFaceColor', cols(1,:))
plot(median(maxRespDiffs(gadNeg)), maxi*1.1, 'v', 'Color', cols(2,:), ...
    'MarkerFaceColor', cols(2,:))
ylim([0 1.2*maxi])
xlabel('Max. response diff.')
ylabel('# Neurons')
xlim(bins([1 end]))
set(gca,'box','off')
legend('Gad+','Gad-')

supprPos = ind & isSuppr==1;
supprNeg = ind & isSuppr==-1;
n1 = hist(maxRespDiffs(supprPos), bins);
n2 = hist(maxRespDiffs(supprNeg), bins);
figure
hold on
plot(bins, n1, 'LineWidth', 2, 'Color', cols(3,:))
plot(bins, n2, 'LineWidth', 2, 'Color', cols(4,:))
maxi = max([n1,n2]);
plot([0 0], [0 maxi.*1.2], 'k--')
plot(median(maxRespDiffs(supprPos)), maxi*1.1, 'v', 'Color', cols(3,:), ...
    'MarkerFaceColor', cols(3,:))
plot(median(maxRespDiffs(supprNeg)), maxi*1.1, 'v', 'Color', cols(4,:), ...
    'MarkerFaceColor', cols(4,:))
ylim([0 1.2*maxi])
xlabel('Max. response diff.')
ylabel('# Neurons')
xlim(bins([1 end]))
set(gca,'box','off')
legend('Suppr+','Suppr-')

%% Compute slopes and intercepts using permutation test
% as far as I remember this didn't produce sensible results
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaFixedAcrossConditions.mat'), 'tuning');
tuning = data.tuning;
numPerms = 500;
fits = struct([]);
for iExp = 1 %:length(tuning)
    fprintf('Dataset %d\n', iExp)
    fits(iExp).subject = tuning(iExp).subject;
    fits(iExp).date = tuning(iExp).date;
    fits(iExp).exp = tuning(iExp).exp;
    fits(iExp).planes = tuning(iExp).planes;
    for iPlane = 1:length(tuning(iExp).plane)
        fprintf('  Plane %d (of %d)\n', iPlane, length(tuning(iExp).plane))
        fprintf('    Cell (of %d):', length(tuning(iExp).plane(iPlane).cellIDs))
        for iCell = 1:length(tuning(iExp).plane(iPlane).cellIDs)
            if mod(iCell,5)==0
                fprintf(' %d', iCell)
            end
            if isempty(tuning(iExp).plane(iPlane).cond(1).cell(iCell).parameters)
                continue
            end
            r1 = tuning(iExp).plane(iPlane).cond(1).cell(iCell).responses;
            r2 = tuning(iExp).plane(iPlane).cond(2).cell(iCell).responses;
            if tuning(iExp).plane(iPlane).isSuppressed(iCell) == 1
                r1 = -r1;
                r2 = -r2;
            end
            r = nansum(cat(3, r1, r2),3);
            resp1 = nanmedian(r1,2);
            resp2 = nanmedian(r2,2);
            maxi = max([resp1; resp2]);
            r = r ./ maxi;
            resp1 = resp1 ./ maxi;
            resp2 = resp2 ./ maxi;
            nonvis = ones(size(r));
            ind = ~isnan(tuning(iExp).plane(iPlane).cond(2).cell(iCell).responses);
            nonvis(ind) = 2;
            numTr = numel(nonvis);
            intercepts = NaN(numPerms,1);
            slopes = NaN(numPerms,1);
            for k = 1:numPerms
                nv = reshape(nonvis(randperm(numTr)),size(nonvis,1), ...
                    size(nonvis,2));
                x1 = NaN(size(r));
                x1(nv==1) = r(nv==1);
                x2 = NaN(size(r));
                x2(nv==2) = r(nv==2);
                med1 = nanmedian(x1,2);
                med2 = nanmedian(x2,2);
                y = fitlm(med1,med2);
                coeffs = y.Coefficients.Estimate;
                intercepts(k) = coeffs(1);
                slopes(k) = coeffs(2);
                
%                 figure
%                 plot(med1,med2,'o')
%                 hold on
%                 plot([0 1],[0 1].*coeffs(2)+coeffs(1),'r')
%                 axis([0 1 0 1])
%                 plot([0 1],[0 1],'k:')
            end
            y = fitlm(resp1,resp2);
            coeffs = y.Coefficients.Estimate;
            fits(iExp).plane(iPlane).cell(iCell).intercept = coeffs(1);
            fits(iExp).plane(iPlane).cell(iCell).slope = coeffs(2);
            fits(iExp).plane(iPlane).cell(iCell).interceptsPermuted = intercepts;
            fits(iExp).plane(iPlane).cell(iCell).slopesPermuted = slopes;
        end
        fprintf('\n')
    end
end
save(fullfile(folderResults, nonVisualSignal, 'slopeFits_lowVsHighNonVisualSignal_permutationTest.mat'),'fits');

%% Plot permutation test results for each cell
data = load(fullfile(folderResults, nonVisualSignal, 'slopeFits_lowVsHighNonVisualSignal_permutationTest.mat'));
fits = data.fits;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaFixedAcrossConditions.mat'), 'tuning');
tuning = data.tuning;
for iExp = 1:length(fits)
    for iPlane = 1:length(fits(iExp).plane)
        fPlots = fullfile(folderResults, ['plots_' nonVisualSignal], ...
            'slopeFits_permutationTest', [tuning(iExp).subject '_' tuning(iExp).date ...
            '_' num2str(tuning(iExp).exp)], num2str(tuning(iExp).planes(iPlane)));
        if ~exist(fPlots, 'dir')
            mkdir(fPlots);
        end
        for iCell = 1:length(fits(iExp).plane(iPlane).cell)
            r1 = tuning(iExp).plane(iPlane).cond(1).cell(iCell).responses;
            r2 = tuning(iExp).plane(iPlane).cond(2).cell(iCell).responses;
            if isempty(r1)
                continue
            end
            if tuning(iExp).plane(iPlane).isSuppressed(iCell) == 1
                r1 = -r1;
                r2 = -r2;
            end
            resp1 = nanmedian(r1,2);
            resp2 = nanmedian(r2,2);
            maxi = max([resp1; resp2]);
            resp1 = resp1 ./ maxi;
            resp2 = resp2 ./ maxi;
            
            figure('Position', [680 680 1045 420])
            subplot(1,2,1)
            plot(resp1, resp2, 'ok')
            hold on
            plot([0 1],[0 1] .* fits(iExp).plane(iPlane).cell(iCell).slope + ...
                fits(iExp).plane(iPlane).cell(iCell).intercept,'r')
            axis([0 1 0 1])
            plot([0 1],[0 1],'k:')
            xlabel(tuning(iExp).plane(iPlane).cond(1).name)
            ylabel(tuning(iExp).plane(iPlane).cond(2).name)
            title(sprintf('Cell %d', iCell))
            
            subplot(1,2,2)
            plot(fits(iExp).plane(iPlane).cell(iCell).interceptsPermuted, ...
                fits(iExp).plane(iPlane).cell(iCell).slopesPermuted, 'k.')
            hold on
            plot(fits(iExp).plane(iPlane).cell(iCell).intercept, ...
                fits(iExp).plane(iPlane).cell(iCell).slope, 'or')
            xlabel('Intercepts')
            ylabel('Slopes')
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(fPlots, sprintf('slopeFit%03d.jpg', iCell)), ...
                '-djpeg','-r0')
            close gcf
        end
    end
end

%% Compute slopes and intercepts using ANOVA
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaFixedAcrossConditions.mat'), 'tuning');
tuning = data.tuning;
fits = struct([]);
fprintf('Dataset: ')
for iExp = 1:length(tuning)
    fprintf('%d ', iExp)
    fits(iExp).subject = tuning(iExp).subject;
    fits(iExp).date = tuning(iExp).date;
    fits(iExp).exp = tuning(iExp).exp;
    fits(iExp).planes = tuning(iExp).planes;
    for iPlane = 1:length(tuning(iExp).plane)
        for iCell = 1:length(tuning(iExp).plane(iPlane).cellIDs)
            if isempty(tuning(iExp).plane(iPlane).cond(1).cell(iCell).parameters)
                continue
            end
            r1 = tuning(iExp).plane(iPlane).cond(1).cell(iCell).responses;
            r2 = tuning(iExp).plane(iPlane).cond(2).cell(iCell).responses;
            r = nansum(cat(3, r1, r2),3);
            stim = repmat((1:size(r,1))',1,size(r,2));
            nonvis = ones(size(r));
            nonvis(~isnan(r2)) = 2;
            pVals = anovan(r(:),[stim(:),nonvis(:)],'model','interaction','display','off');
            fits(iExp).plane(iPlane).cell(iCell).ANOVApStim = pVals(1);
            fits(iExp).plane(iPlane).cell(iCell).ANOVApNonvis = pVals(2);
            fits(iExp).plane(iPlane).cell(iCell).ANOVApInteract = pVals(3);
            
            resp1 = nanmedian(r1,2);
            resp2 = nanmedian(r2,2);
            if tuning(iExp).plane(iPlane).isSuppressed(iCell) == 1
                resp1 = -resp1;
                resp2 = -resp2;
            end
            maxi = max([resp1; resp2]);
            resp1 = resp1 ./ maxi;
            resp2 = resp2 ./ maxi;
            if all(pVals >= 0.05) % neuron not modulated by stimulus or nonvisual signal
                fits(iExp).plane(iPlane).cell(iCell).adjR2 = [];
                fits(iExp).plane(iPlane).cell(iCell).intercept = [];
                fits(iExp).plane(iPlane).cell(iCell).interceptCI = [];
                fits(iExp).plane(iPlane).cell(iCell).slope = [];
                fits(iExp).plane(iPlane).cell(iCell).slopeCI = [];
                fits(iExp).plane(iPlane).cell(iCell).maxRespDiff = [];
                fits(iExp).plane(iPlane).cell(iCell).modulations = 'none';
            elseif pVals(2) < 0.05 && pVals(3) >= 0.05 % purely additive
                [y,gof] = fit(resp1, resp2, @(a,x) x+a, ...
                    'StartPoint', median(resp1)-median(resp2));
                fits(iExp).plane(iPlane).cell(iCell).adjR2 = gof.adjrsquare;
                fits(iExp).plane(iPlane).cell(iCell).intercept = y.a;
                fits(iExp).plane(iPlane).cell(iCell).interceptCI = confint(y)';
                fits(iExp).plane(iPlane).cell(iCell).slope = 1;
                fits(iExp).plane(iPlane).cell(iCell).slopeCI = [];
                fits(iExp).plane(iPlane).cell(iCell).maxRespDiff = y.a;
                fits(iExp).plane(iPlane).cell(iCell).modulations = 'additive';
%             elseif pVals(2) >= 0.05 && pVals(3) < 0.05 % purely multiplicative
%                 meanResp = nanmean([resp1; resp2]);
%                 resp1 = resp1 - meanResp; % subtract mean response so that no intercept is necessary when fitting
%                 resp2 = resp2 - meanResp;
%                 [y,gof] = fit(resp1, resp2, @(a,x) x.*a, ...
%                     'StartPoint', nanmean(resp1./resp2));
%                 intercept = meanResp - y.a * meanResp; % shift fitted line by meanResp on x and y axis
%                 % calculate response difference;
%                 fits(iExp).plane(iPlane).cell(iCell).adjR2 = gof.adjrsquare;
%                 fits(iExp).plane(iPlane).cell(iCell).intercept = intercept;
%                 fits(iExp).plane(iPlane).cell(iCell).interceptCI = [];
%                 fits(iExp).plane(iPlane).cell(iCell).slope = y.a;
%                 fits(iExp).plane(iPlane).cell(iCell).slopeCI = confint(y)';
%                 % calculate response difference
%                 if abs(y.a) <= 1
%                     % if the absolute slope is smaller one, then the response
%                     % to the most driving stimulus during small pupil is 1
%                     respDiff = y.a+intercept - 1;
%                 else
%                     % if the absolute slope is larger one, then the response
%                     % to the most driving stimulus during large pupil is 1
%                     respDiff = 1 - (1-intercept)/y.a;
%                 end
%                 fits(iExp).plane(iPlane).cell(iCell).maxRespDiff = respDiff;
%                 if gof.adjrsquare < 0 % check whether the interaction is multiplicative
%                     fits(iExp).plane(iPlane).cell(iCell).modulations = 'none';
%                 else
%                     fits(iExp).plane(iPlane).cell(iCell).modulations = 'multiplicative';
%                 end
            elseif pVals(2) < 0.05 || pVals(3) < 0.05 % not purely additive
                y = fitlm(resp1,resp2);
                coefficients = y.Coefficients.Estimate;
                CI = y.coefCI;
                fits(iExp).plane(iPlane).cell(iCell).adjR2 = y.Rsquared.Adjusted;
                fits(iExp).plane(iPlane).cell(iCell).intercept = coefficients(1);
                fits(iExp).plane(iPlane).cell(iCell).interceptCI = CI(1,:);
                fits(iExp).plane(iPlane).cell(iCell).slope = coefficients(2);
                fits(iExp).plane(iPlane).cell(iCell).slopeCI = CI(2,:);
                if abs(coefficients(2)) <= 1
                    respDiff = coefficients(2)+coefficients(1) - 1;
                else
                    respDiff = 1 - (1-coefficients(1))/coefficients(2);
                end
                fits(iExp).plane(iPlane).cell(iCell).maxRespDiff = respDiff;
                if y.Rsquared.Adjusted < 0
                    fits(iExp).plane(iPlane).cell(iCell).modulations = 'none';
                else
                    fits(iExp).plane(iPlane).cell(iCell).modulations = 'mixed';
                end
            else % only modulated by stimulus 
                fits(iExp).plane(iPlane).cell(iCell).adjR2 = [];
                fits(iExp).plane(iPlane).cell(iCell).intercept = [];
                fits(iExp).plane(iPlane).cell(iCell).interceptCI = [];
                fits(iExp).plane(iPlane).cell(iCell).slope = [];
                fits(iExp).plane(iPlane).cell(iCell).slopeCI = [];
                fits(iExp).plane(iPlane).cell(iCell).maxRespDiff = [];
                fits(iExp).plane(iPlane).cell(iCell).modulations = 'stimOnly';
            end
        end
    end
end
fprintf('\n')
save(fullfile(folderResults, nonVisualSignal, ...
    'slopeFits_lowVsHighNonVisualSignal.mat'),'fits');

%% Plot responses differences against stimulus modulation and pref. ori.
% load fits (containing response differences)
data = load(fullfile(folderResults, nonVisualSignal, ...
    'slopeFits_lowVsHighNonVisualSignal_onlySignifParamsFit.mat'));
fits = data.fits;
% load tuning (containing pref. oris)
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaFixedAcrossConditions_simultaneous.mat'));
tuning = data.tuning;
% load results (containing explained variance by stimulus alone)
data = load(fullfile(folderResults, 'linFit_running_pupil_explainedVars.mat'));
results = data.results;

respDiffs = [];
explVarStim = [];
prefDir = [];
isSuppr = [];
isGad = [];
stimModel = 12;
dirParam = 1;
for iExp = 1:length(fits)
    for iPlane = 1:length(fits(iExp).plane)
        explVars = {results(iExp).plane(iPlane).crossVal.explVars};
        ev = NaN(length(explVars),1);
        ind = ~cellfun(@isempty, explVars);
        explVars = cat(1,explVars{:});
        ev(ind) = explVars(:,12);
        explVarStim = [explVarStim; ev];
        
        rDiffs = {fits(iExp).plane(iPlane).cell.maxRespDiff};
        rd = NaN(size(ev));
        ind = ~cellfun(@isempty, rDiffs);
        rd(ind) = cell2mat(rDiffs);
        respDiffs = [respDiffs; rd];
        
        params = {tuning(iExp).plane(iPlane).cond(1).cell.parameters};
        pars = NaN(size(ev));
        ind = ~cellfun(@isempty, params);
        params = cell2mat(params);
        pars(ind) = params(dirParam,:);
        prefDir = [prefDir; pars];
        
        isSuppr = [isSuppr; tuning(iExp).plane(iPlane).isSuppressed];
        isGad = [isGad; tuning(iExp).plane(iPlane).isGad];
    end
end

figure
plot(explVarStim, abs(respDiffs), 'k.')
ind = all(~isnan([explVarStim, respDiffs]),2);
[rho,p] = corr(explVarStim(ind), abs(respDiffs(ind)));
xlabel('Expl. var. by stimulus')
ylabel('Response difference (large vs small pupil)')
title(sprintf('Rho = %.3f (p = %.3f)',rho,p))
set(gca,'box','off')

figure
plot(prefDir,abs(respDiffs), 'k.')
set(gca,'box','off','XTick',0:180:360)
xlabel('Direction')
ylabel('Response difference (large vs small pupil)')
xlim([0 360])

%% Compare tuning curves across conditions
% Compare tuning properties across nonvisual conditions
R2 = [];
isSuppressed = [];
isGadAdditive = [];
for iExp = 1:length(tuning)
    for iPlane = 1:length(tuning(iExp).plane)
        R2 = [R2; tuning(iExp).plane(iPlane).R2];
        isSuppressed = [isSuppressed; tuning(iExp).plane(iPlane).isSuppressed];
        isGadAdditive = [isGadAdditive; tuning(iExp).plane(iPlane).isGad];
    end
end
ind = isnan(R2);
R2(ind) = [];
isSuppressed(ind) = [];
isGadAdditive(ind) = [];

% response at pref. dir. (defined during low behav. value!!!)
maxResp = [];
minResp = [];
respAtPref = [];
directionIndex = [];
fullWidthHalfHeight = [];
blankResp = [];
medianResp = [];
for iExp = 1:length(tuning)
    for iPlane = 1:length(tuning(iExp).plane)
        pars1 = cat(2,tuning(iExp).plane(iPlane).cond(1).cell.parameters);
        pars2 = cat(2,tuning(iExp).plane(iPlane).cond(2).cell.parameters);
        
        curves = cell(1,2);
        curves{1} = cat(1,tuning(iExp).plane(iPlane).cond(1).cell.curve);
        curves{2} = cat(1,tuning(iExp).plane(iPlane).cond(2).cell.curve);
        
        maxResp = [maxResp; ...
            [sum(pars1([2,4],:),1)', sum(pars2([2,4],:),1)']];
        minR1 = min(curves{1},[],2);
        minEV = min(curves{2},[],2);
        minResp = [minResp; [minR1, minEV]];
        Rn1 = (1-pars1(3,:))./(1+pars1(3,:)) ./ pars1(2,:);
        Rn2 = (1-pars2(3,:))./(1+pars2(3,:)) ./ pars2(2,:);
        r = [sum(pars1([2,4],:),1)', sum(pars2([2,4],:),1)'];
        ind = abs(pars1(1,:) - pars2(1,:)) > 90;
        r(ind,2) = Rn2(ind) + pars2(4,ind);
        respAtPref = [respAtPref; r];
        r = [((pars1(2,:)-Rn1)./(pars1(2,:)+Rn1-2*minR1'))', ...
            ((pars2(2,:)-Rn2)./(pars2(2,:)+Rn2-2*minEV'))'];
        directionIndex = [directionIndex; r];
        
        for j = 1:size(curves{1},1)
            fwhh = zeros(1,2);
            for c = 1:2
                [maxi,indMax] = max(curves{c}(j,:));
                [mini,indMin] = min(curves{c}(j,:));
                if maxi == mini
                    fwhh(c) = NaN;
                    continue
                end
                halfHeight = mini+.5*(maxi-mini);
                curv = [curves{c}(j,indMax:end),curves{c}(j,:)];
                ind = find(curv < halfHeight,1);
                if abs(curv(ind-1)-halfHeight)<abs(curv(ind)-halfHeight)
                    ind = ind-1;
                end
                fwhh(c) = (2*ind-1) * diff(x([1 2]));
            end
            fullWidthHalfHeight = [fullWidthHalfHeight; fwhh];
        end
        
        blankResp = [blankResp; [mean(cat(1,tuning(iExp).plane(iPlane) ...
            .cond(1).cell.blankResponses),2), ...
            mean(cat(1,tuning(iExp).plane(iPlane) ...
            .cond(2).cell.blankResponses),2)]];
        med1 = cat(3,tuning(iExp).plane(iPlane).cond(1).cell.responses);
        med2 = cat(3,tuning(iExp).plane(iPlane).cond(2).cell.responses);
        medianResp = [medianResp; [squeeze(nanmedian(reshape( ...
            med1,[],1,size(med1,3)),1)), ...
            squeeze(nanmedian(reshape(med2,[],1,size(med2,3)),1))]];
    end
end

isBehaviorEnhanced = double(maxResp(:,1) < maxResp(:,2));
measures = {maxResp, minResp, respAtPref, directionIndex, ...
    fullWidthHalfHeight, blankResp};
titles = {'Maximal response','Minimal response', ...
    'Response at pref. direction (default condition)', ...
    'Direction index', ...
    'Full width at half height (around each pref. dir.)', ...
    'Response to blank screen'};
conditions = {isSuppressed, isBehaviorEnhanced, isGadAdditive};
condValues = [-1 1; 1 0; -1 1];
condLabels = {'Contrast enhanced','Contrast suppressed'; ...
    'Behavior enhanced','Behavior suppressed'; ...
    'Excitatory','Inhibitory'};

% plot maximum responses (high-low beh. value) vs. minimum response
% (high-low beh. value)
ind = R2>=.3;
figure
mx = maxResp(ind,2)-maxResp(ind,1);
my = minResp(ind,2)-minResp(ind,1);
scatter(mx,my)
hold on
mini = min([mx;my]);
maxi = max([mx;my]);
plot([mini maxi],[mini maxi],'r')
plot(median(mx),median(my),'+r','LineWidth',2)
axis([mini maxi mini maxi])
axis square
title('Change in max vs. min response')
xlabel(sprintf('Max (%s - %s)',labels{2},labels{1}))
ylabel(sprintf('Min (%s - %s)',labels{2},labels{1}))

% plot maximum responses (high-low beh. value) vs. minimum response
% (high-low beh. value) separately for each cell group
figure('Position',[1929 2 756 1114])
for cond = 1:length(conditions)
    for v = 1:2
        cnd = conditions{cond}(R2>=0.3);
        subplot(3,2,(cond-1)*2+v)
        scatter(mx(cnd==condValues(cond,v)), ...
            my(cnd==condValues(cond,v)))
        hold on
        mini = min([mx(cnd==condValues(cond,v));my(cnd==condValues(cond,v))]);
        maxi = max([mx(cnd==condValues(cond,v));my(cnd==condValues(cond,v))]);
        plot([mini maxi],[mini maxi],'r')
        axis([mini maxi mini maxi])
        axis square
        title(condLabels{cond,v})
        xlabel(sprintf('Max (%s - %s)',labels{2},labels{1}))
        ylabel(sprintf('Min (%s - %s)',labels{2},labels{1}))
    end
end
annotation('textbox',[.3 .95 .4 .04],'String','Change in max vs. min response', ...
    'FontSize',12,'LineStyle','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle')

for m = 1:length(measures)
    meas = measures{m};
    if isempty(meas)
        continue
    end
    meas(R2<0.3,:) = [];
    figure
    scatter(meas(:,1), meas(:,2))
    hold on
    mini = min(meas(:));
    maxi = max(meas(:));
    plot([mini maxi],[mini maxi],'r')
    axis([mini maxi mini maxi])
    axis square
    title(titles{m})
    xlabel(labels{1})
    ylabel(labels{2})
    
    figure('Position',[1929 2 756 1114])
    for cond = 1:length(conditions)
        for v = 1:2
            cnd = conditions{cond}(R2>=0.3);
            subplot(3,2,(cond-1)*2+v)
            scatter(meas(cnd==condValues(cond,v),1), ...
                meas(cnd==condValues(cond,v),2))
            hold on
            mini = min(reshape(meas(cnd==condValues(cond,v),:),[],1));
            maxi = max(reshape(meas(cnd==condValues(cond,v),:),[],1));
            plot([mini maxi],[mini maxi],'r')
            plot(nanmedian(meas(cnd==condValues(cond,v),1)), ...
                nanmedian(meas(cnd==condValues(cond,v),2)), 'r+', 'LineWidth',2)
            axis([mini maxi mini maxi])
            axis square
            title(condLabels{cond,v})
            xlabel(labels{1})
            ylabel(labels{2})
        end
    end
    annotation('textbox',[.3 .95 .4 .04],'String',titles{m}, ...
        'FontSize',12,'LineStyle','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle')
end

meas = medianResp(R2<.3,:);
figure
scatter(meas(:,1),meas(:,2))
hold on
mini = min(reshape(meas,[],1));
maxi = max(reshape(meas,[],1));
plot([mini maxi],[mini maxi],'r')
axis([mini maxi mini maxi])
axis square
title('Median visual response of untuned neurons')
xlabel(labels{1})
ylabel(labels{2})
figure('Position',[1929 300 756 809])
conds = [1 3];
for c = 1:2
    cond = conds(c);
    for v = 1:2
        cnd = conditions{cond}(R2<0.3);
        subplot(2,2,(c-1)*2+v)
        scatter(meas(cnd==condValues(cond,v),1), ...
            meas(cnd==condValues(cond,v),2))
        hold on
        mini = min(reshape(meas(cnd==condValues(cond,v),:),[],1));
        maxi = max(reshape(meas(cnd==condValues(cond,v),:),[],1));
        plot([mini maxi],[mini maxi],'r')
        axis([mini maxi mini maxi])
        axis square
        title(condLabels{cond,v})
        xlabel(labels{1})
        ylabel(labels{2})
    end
end
annotation('textbox',[.3 .95 .4 .04],'String', ...
    'Median visual response of untuned neurons', ...
    'FontSize',12,'LineStyle','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle')
