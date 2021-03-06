%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish';
folderTools = 'C:\STORAGE\workspaces';
folderThisRepo = 'C:\dev\workspace\schroeder-et-al-2020';

%% Parameters
% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

% to group into ON, OFF and ON+OFF
onThr = 0.5; % ON field is at least three times stronger than OFF field
offThr = -0.5; % OFF field is at least three times stronger than ON field

% colormaps
red = [1 0 .5];
blue = [0 .5 1];
black = [1 1 1].*0.5;
grad = linspace(0,1,100)';
reds = red.*flip(grad) + [1 1 1].*grad;
blacks = black.*flip(grad) + [1 1 1].*grad;
cm_ON = [blacks; flip(reds(1:end-1,:),1)];
blues = blue.*flip(grad) + [1 1 1].*grad;
cm_OFF = [blacks; flip(blues(1:end-1,:),1)];
colormaps = {cm_ON, cm_OFF};

%% Examples
examplesRF = {'SS078', '2017-10-05', 1, 156, -1; ...
            'SS078', '2017-10-05', 1, 225, 1; ...
            'SS078', '2017-10-05', 1, 22, -1; ...
            'SS078', '2017-09-28', 1, 243, 1; ...
            'SS069', '2016-10-21', 1, 223, -1; ...
            'SS078', '2017-09-28', 1, 263, 1};

examplesTun = {'SS076', '2017-10-04', 1, 137; ...
            'SS077', '2017-10-05', 1, 24; ...
            'SS069', '2016-10-13', 3, 4};

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Figure 1F (receptive fields of example boutons)
titles = {'ON field','OFF field'};
for ex = 1:size(examplesRF,1)
    roiPlanes = readNPY(fullfile(folderBase, 'boutons', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_2pRois._ss_2pPlanes.npy'));
    ids = readNPY(fullfile(folderBase, 'boutons', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_2pRois.ids.npy'));
    rfs = readNPY(fullfile(folderBase, 'boutons', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_rf.maps.npy'));
    pos = readNPY(fullfile(folderBase, 'boutons', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_rfDescr.edges.npy'));
    
    roi = roiPlanes == examplesRF{ex,3} & ids == examplesRF{ex,4};
    rf = squeeze(rfs(roi,:,:,:,:));
    rf(:,:,:,2) = -rf(:,:,:,2);
    squW = diff(pos(1:2)) / size(rf,2);
    squH = diff(pos(3:4)) / size(rf,1);
    [mx,mxTime] = max(max(reshape(permute(abs(rf),[1 2 4 3]),[],size(rf,3)),[],1));
    % fit Gaussian
    fitPars = whiteNoise.fit2dGaussRF(mean(rf(:,:,mxTime,:),4), false);
    [x0, y0] = meshgrid(linspace(0.5, size(rf,2)+0.5, 100), ...
        linspace(0.5, size(rf,1)+0.5, 100));
    fitRF = whiteNoise.D2GaussFunctionRot(fitPars, cat(3, x0, y0));
    [x1,y1] = meshgrid(linspace(pos(1), pos(2), 100), ...
        linspace(pos(3), pos(4), 100));
    figure
    for f = 1:2
        subplot(2,1,f)
        imagesc([pos(1)+squW/2 pos(2)-squW/2], [pos(3)+squH/2 pos(4)-squH/2], ...
            rf(:,:,mxTime,f),[-mx mx])
        hold on
        contour(x1, y1, fitRF, [1 1].*fitPars(1)/2, 'k', 'LineWidth', 1)
        axis image
        set(gca, 'box', 'off', 'XTick', pos(1:2), 'YTick', [pos(3) 0 pos(4)], ...
            'YTickLabel', [-pos(3) 0 -pos(4)])
        colormap(gca, colormaps{f})
        ylabel(titles{f})
        colorbar
        if pos(2) > 0 % RF map crosses vertical midline
            xlim([pos(1) 0])
            set(gca, 'XTick', [pos(1) 0])
        end
        if f == 1
            title(sprintf('%s, %s, plane %d, unit %d', ...
                examplesRF{ex,1}, examplesRF{ex,2}, examplesRF{ex,3}, examplesRF{ex,4}), ...
                'interpreter', 'none')
        end
    end
end

%% Figure 1G (histogram of ON/OFF indices)
% Collect relevant variables from visual noise data
subjects = {};
dates = {};
planes = [];
ids = [];
evStim = [];
lambdasStim = [];
pValues = [];
OnOffValues = [];

subjDirs = dir(fullfile(folderBase, 'boutons', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'boutons', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_rf.maps.npy'))
            continue
        end
        
        planes = [planes; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_2pRois._ss_2pPlanes.npy'))];
        ids = [ids; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_2pRois.ids.npy'))];
        rfs_exp = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_rf.maps.npy'));
        evStim = [evStim; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_rf.explVarsStim.npy'))];
        lambdasStim = [lambdasStim; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_rf.lambdasStim.npy'))];
        pValues = [pValues; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_rf.pValues.npy'))];
        
        n = size(rfs_exp,1);
        subjects = [subjects; repmat({name}, n, 1)];
        dates = [dates; repmat({date}, n, 1)];
        
        oov = NaN(n,2);
        for iCell = 1:n
            rf = squeeze(rfs_exp(iCell,:,:,:,:));
            [~,t] = max(max(reshape(permute(abs(rf),[1 2 4 3]), [], ...
                size(rf,3)), [], 1));
            rf = squeeze(rf(:,:,t,:));
            [mx,row] = max(max(abs(rf),[],3),[],1);
            [~,col] = max(mx);
            row = row(col);
            rf = squeeze(rf(row,col,:));
            rf(2) = -rf(2);
            oov(iCell,:) = rf;
        end
        OnOffValues = [OnOffValues; oov];
    end
end

[~,type] = max(abs(OnOffValues),[],2);
ind = sub2ind(size(OnOffValues), (1:size(OnOffValues,1))', type);
signs = sign(OnOffValues(ind));
OnOffRatios = OnOffValues;
OnOffRatios(signs<0,:) = -OnOffRatios(signs<0,:);
OnOffRatios(OnOffRatios<0) = 0;
OnOffRatios = (OnOffRatios(:,1)-OnOffRatios(:,2))./sum(OnOffRatios,2);

validRF = pValues < minPVal & evStim > minExplainedVarianceStim & ...
    lambdasStim < maxLambda;

% Plot distribution of ON-OFF-ratios
binSize = 0.05;
edges = -1:binSize:1;
bins = edges(1:end-1)+binSize/2;
figure
histogram(OnOffRatios(validRF), edges, 'FaceColor', 'k')
n = histcounts(OnOffRatios(validRF), edges);
mx = max(n) * 1.05;
hold on
plot([1 1].*offThr, [0 mx], 'k')
plot([1 1].*onThr, [0 mx], 'k')
for ex = 1:size(examplesRF,1)
    idx = strcmp(subjects, examplesRF{ex,1}) & strcmp(dates, examplesRF{ex,2}) & ...
        planes==examplesRF{ex,3} & ids==examplesRF{ex,4};
    if OnOffRatios(idx) > onThr
        color = red;
    elseif OnOffRatios(idx) < offThr
        color = blue;
    else
        color = (red+blue)./2;
    end
    if examplesRF{ex,5} > 0
        color = color .* 0.7 + [1 1 1].*0.3;
    else
        color = color .* 0.7 + [0 0 0].*0.3;
    end
    
    plot(OnOffRatios(idx), mx, 'v', ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none')        
end
set(gca, 'box', 'off', 'XTick', [-1 offThr 0 onThr 1])
xlim([-1 1])
ylim([0 mx*1.05])
xlabel('ON/OFF index')
ylabel('#Boutons')
title(sprintf('n = %d', sum(validRF)))

%% Figure 1H & I (responses to single gratings and tuning curves for example boutons)
buffer = 2; % in sec (before and after stim period)
for ex = 1:size(examplesTun,1)
    planes = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pRois._ss_2pPlanes.npy'));
    ids = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pRois.ids.npy'));
    planeDelays = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pPlanes.delay.npy'));
    traces = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pCalcium.dff.npy'));
    time = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pCalcium.timestamps.npy'));
    predictions = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingPredictions.dff.npy'));
    predTime = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingPredictions.timestamps.npy'));
    stimIntervals = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_grating.intervals.npy'));
    stimSequence = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_grating._ss_gratingID.npy'));
    directions = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingID.directions.npy'));
    largePupil = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingTrials.largePupil.npy'));
    curves = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_tuning.curvesSmall.npy'));
    amplitudes = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingTrials.amplitudes.npy'));
    
    cellID = examplesTun{ex,3}==planes & examplesTun{ex,4}==ids;
    t = time + planeDelays(examplesTun{ex,3});
    stimMatrix = exp.buildStimMatrix(stimSequence, stimIntervals, t);
    directions(isnan(directions)) = [];
    timeBin = median(diff(time));
    repetitions = sum(stimSequence == 1);
    stimDurInFrames = round(sum(stimMatrix(1,:)) / repetitions);
    stimDur = stimDurInFrames * timeBin;
    offset = ceil(buffer / timeBin);
    resp = squeeze(exp.getTracesPerStimulus(traces(:,cellID), stimMatrix, ...
        [1 1] .* offset)); % [stimulus x trial x time]
    predResp = squeeze(exp.getTracesPerStimulus(predictions(:,cellID), ...
        stimMatrix, [1 1] .* offset));
    
    ind = ~largePupil;
    ind = repmat(ind, 1, 1, size(resp,3));
    temp = resp(1:end-1,:,:);
    temp(~ind) = NaN;
    respMean = nanmean(temp, 2);
    temp = predResp(1:end-1,:,:);
    temp(~ind) = NaN;
    predMean = nanmean(temp, 2);
    respMean = respMean(1:3:end,:);
    predMean = predMean(1:3:end,:);
    
    mini = min([respMean(:); predMean(:)]);
    maxi = max([respMean(:); predMean(:)]);
    rng = maxi - mini;
    mini = mini - 0.05*rng;
    maxi = maxi + 0.05*rng;
    xDist = .5;
    traceDur = stimDur + 2*buffer;
    respTime = (-offset:stimDurInFrames+offset-1) .* timeBin;
    
    figure('Position',[50 500 1000 420])
    hold on
    h = [0 0];
    x0 = 0;
    for st = 1:size(respMean,1)
        fill([0 stimDur stimDur 0] + x0, ...
            [mini mini maxi maxi], 'k', 'FaceColor', 'k', ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none')
        plot(respTime([1 end]) + x0, [0 0], 'k:')
        h1 = plot(respTime + x0, squeeze(respMean(st,:)), ...
            'Color', 'k');
        h2 = plot(respTime + x0, squeeze(predMean(st,:)), ...
            'Color', 'k', 'LineWidth', 2);
        x0 = x0 + traceDur + xDist;
        if st==1
            h(1) = h2;
            h(2) = h1;
        end
    end
    axis tight
    set(gca, 'XTick', [0 stimDur])
    legend(h, {'Prediction from kernel fit','Data'})
    xlabel('Stimuli')
    ylabel('\DeltaF/F')
    title(sprintf('Example %d',ex))
    
    figure
    hold on
    amps = amplitudes(:,:,cellID);
    amps(largePupil) = NaN;
    plot(1:360, curves(cellID,:), 'Color', 'k', 'LineWidth',2);
    m = nanmean(amps,2);
    s = nanstd(amps,0,2) ./ sqrt(sum(~largePupil,2));
    errorbar([directions; 360], m([1:end 1]), s([1:end 1]), 'o', ...
        'Color', 'k', 'CapSize', 2, 'MarkerFaceColor', 'k')
    plot([0 360], [0 0], 'k:', 'LineWidth', 2)
    set(gca,'XTick',0:90:360)
    title(sprintf('Example %d',ex))
    xlim([0 360])
    mini = min([mini; m-s; curves(cellID,:)']);
    maxi = max([maxi; m+s; curves(cellID,:)']);
    ylim([mini maxi])
    xlabel('Direction (in degrees)')
    ylabel('\DeltaF/F')
end

%% Prepare for Figs. 1J-L
% Collect relevant variables from tuning data
subjects = {};
dates = {};
planes = [];
ids = [];
prefDirs = [];
minima = [];
maxima = [];
means = [];
isTuned = [];
isSuppr = [];
amplitudes = {};
largePupil = {};

subjDirs = dir(fullfile(folderBase, 'boutons', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'boutons', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.parametersSmall.npy'))
            continue
        end
        
        p = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_2pRois._ss_2pPlanes.npy'));
        n = length(p);
        planes = [planes; p];
        ids = [ids; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_2pRois.ids.npy'))];
        subjects = [subjects; repmat({name}, n, 1)];
        dates = [dates; repmat({date}, n, 1)];
        parsS = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.parametersSmall.npy'));
        parsL = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.parametersLarge.npy'));
        isT = [~isnan(parsS(:,2)), ~isnan(parsL(:,2))];
        isTuned = [isTuned; all(isT,2)];
        prefD = parsS(:,1);
        prefD(~isT(:,1)) = NaN;
        prefDirs = [prefDirs; prefD];
        curves = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.curvesSmall.npy'));
        isSuppr = [isSuppr; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.isSuppressed.npy'))];
        amp = readNPY(fullfile(folderBase, ...
            'boutons', name, date, '001\_ss_gratingTrials.amplitudes.npy'));
        amplitudes = [amplitudes; permute(mat2cell(amp, ...
            size(amp,1), size(amp,2), ones(1,n)), [3 1 2])];
        largePupil = [largePupil; repmat({readNPY(fullfile(folderBase, ...
            'boutons', name, date, '001\_ss_gratingTrials.largePupil.npy'))}, n, 1)];
        
        mi = NaN(n,1);
        ma = NaN(n,1);
        mn = NaN(n,1);
        for iCell = 1:n
            if isT(iCell,1) % tuned when pupil small
                mn(iCell) = mean(curves(iCell,:));
                oris = mod(prefD(iCell) + [0 90 180], 360);
                respS = gratings.orituneWrappedConditions(parsS(iCell,:), oris);
                ma(iCell) = respS(1);
                if respS(1)-respS(2) < 0 % suppressed by gratings
                    mi(iCell) = max(respS(2:3));
                else
                    mi(iCell) = min(respS(2:3));
                end
            else
                mn(iCell) = parsS(iCell,1);
            end
        end
        minima = [minima; mi];
        maxima = [maxima; ma];
        means = [means; mn];
    end
end

shuffles = 1000;
directions = 0:30:330;
dirVectors = exp(directions./180.*pi .* 1i);
oriVectors = exp(directions./180.*2.*pi .* 1i);

ampSmall = amplitudes; % {nROIs x 1}, each entry: [nStim x nReps]
ampShuffle = cell(size(amplitudes));  % {nROIs x 1}, each entry: [nReps x nStim x nShuffles]
sz = cell2mat(cellfun(@size, amplitudes, 'UniformOutput', false));
numReps = unique(sz(:,2));
numStim = unique(sz(:,1));
permutations = cell(1, length(numReps));
for m = 1:length(numReps)
    permutations{m} = NaN(numReps(m)*numStim,1);
    for sh = 1:shuffles
        permutations{m}(:,sh) = randperm(numReps(m)*numStim);
    end
end
for n = 1:size(amplitudes,1)
    m = find(numReps == size(amplitudes{n},2));
    if all(isnan(ampSmall{n}(:)))
        ampShuffle{n} = NaN(numReps(m), numStim, shuffles);
        continue
    end
    ampSmall{n}(largePupil{n}) = NaN;
    a = reshape(ampSmall{n}, [], 1);
    a = a(permutations{m});
    ampShuffle{n} = reshape(a, numReps(m), numStim, shuffles);
end
meanAmp = cellfun(@nanmean, ampSmall, repmat({2},size(amplitudes,1),1), ...
    'UniformOutput', false);
meanAmp = cell2mat(meanAmp')'; % [ROIs x stimuli], amplitudes averaged across small pupil trials
shuffleAmp = cellfun(@nanmean, ampShuffle, 'UniformOutput', false);
shuffleAmp = cell2mat(shuffleAmp); % [ROIs x stimuli x shuffles] mean amplitudes after shuffling stimulus labels

% inverte responses of suppressed ROIs
meanAmp(isSuppr==1,:) = -meanAmp(isSuppr==1,:);
shuffleAmp(isSuppr==1,:,:) = -shuffleAmp(isSuppr==1,:,:);
% set responses below baseline to zero
meanAmp(meanAmp<0) = 0;
shuffleAmp(shuffleAmp<0) = 0;
% standardize responses so they sum to one
meanAmp = meanAmp ./ nansum(meanAmp,2);
shuffleAmp = shuffleAmp ./ nansum(shuffleAmp,2);

% Determine DSIs
vects = sum(dirVectors .* meanAmp, 2);
shuffleVects = squeeze(sum(dirVectors .* shuffleAmp, 2));
DSIs = abs(vects);
nullDSIs = abs(shuffleVects);
p_DSI = sum(nullDSIs > DSIs,2) ./ shuffles;

% Determine OSIs
vects = sum(oriVectors .* meanAmp, 2);
shuffleVects = squeeze(sum(oriVectors .* shuffleAmp, 2));
OSIs = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSI = sum(nullOSIs > OSIs,2) ./ shuffles;

%% Figure 1J (histogram of maximum amplitudes in response to gratings)
cols = lines(3);
xLim = [9 10];
yLim = [3 5];

mi = minima;
ma = maxima;
j = isnan(mi); % not tuned when pupil small
mi(j) = means(j);
ma(j) = means(j);

j = ma - mi < 0; % if suppressed by gratings, swap responses to preferred
                 % and non-preferred stimulus because response to pref is
                 % minimum response
m = mi(j);
mi(j) = ma(j);
ma(j) = m;
b_mi = [-flip(2.^(-3:xLim(1))), 0, 2.^(-3:xLim(2))];
b_ma = [-flip(2.^(-3:yLim(1))), 0, 2.^(-3:yLim(2))];

figure
N1 = histcounts(ma(isTuned==1), b_ma);
N2 = histcounts(ma(isTuned==0), b_ma);
b = bar(1.5:length(b_ma), [N1',N2'], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
hold on
for ex = 1:size(examplesTun,1)
    idx = strcmp(subjects, examplesTun{ex,1}) & ...
        strcmp(dates, examplesTun{ex,2}) & planes==examplesTun{ex,3} & ...
        ids==examplesTun{ex,4};
    exMax = find(b_ma<ma(idx),1,'last');
    if ma(idx)>0
        exMax = exMax + rem(log2(ma(idx)),1);
    else
        exMax = exMax - rem(log2(-ma(idx)),1);
    end
    plot(exMax,520,'v','Color',cols(ex,:))
end
xlabel('Maximum')
ylabel('#Boutons')
title(sprintf('n = %d', sum(~isnan(ma))))
set(gca,'XTick',1:length(b_ma),'XTickLabel',b_ma,'box','off')
legend(b,{'tuned','not tuned'},'Location','NorthWest')

%% Figure 1K (scatter of DSIs vs OSIs)
% Test significance for separation between high OSIs and high DSIs
goodDSIs = DSIs(~isnan(DSIs) & isTuned & (p_DSI<.05|p_OSI<.05));
goodOSIs = OSIs(~isnan(OSIs) & isTuned & (p_DSI<.05|p_OSI<.05));
numPerm = 10000;
permutations = NaN(length(goodDSIs), numPerm);
for p = 1:size(permutations,2)
    permutations(:,p) = randperm(size(permutations,1));
end
permOSIs = goodOSIs(permutations);
meanAbsDiff = mean(abs(goodDSIs - goodOSIs));
permAbsDiffs = mean(abs(goodDSIs - permOSIs), 1);
pDiff = sum(permAbsDiffs > meanAbsDiff) / numPerm;

% Plot scatter of DSI versus OSI
cols = lines(4);
binSize = 0.02;
mx = ceil(max([DSIs; OSIs]) / binSize) * binSize;
colors = [1 0 0; 0 0 1; .5 0 .5];
inds = [p_DSI<.05&p_OSI>=.05, p_DSI>=.05&p_OSI<.05, p_DSI<.05&p_OSI<.05];
h = [0 0 0];
figure
hold on
plot([0 mx], [0 mx], 'k')
for k = 3:-1:1
    h(k) = scatter(DSIs(isTuned & inds(:,k)), ...
        OSIs(isTuned & inds(:,k)), [], colors(k,:), 'filled');
end
for ex = 1:size(examplesTun,1)
    idx = strcmp(subjects, examplesTun{ex,1}) & ...
        strcmp(dates, examplesTun{ex,2}) & planes==examplesTun{ex,3} & ...
        ids==examplesTun{ex,4};
    if p_DSI(idx)>=0.05 && p_OSI(idx)>=0.05
        continue
    end
    scatter(DSIs(idx), OSIs(idx), [], cols(ex,:), 'LineWidth', 2)
end
legend(h,'dir-sel','ori-sel','dir&ori')
axis square
axis([0 mx 0 mx])
set(gca, 'box','off','XTick',0:.1:mx,'YTick',0:.1:mx)
xlabel('DSI')
ylabel('OSI')
title(sprintf('n = %d (DSI and OSI are different, p = %.2e)', ...
    sum(~isnan(DSIs) & isTuned & (p_DSI<0.05|p_OSI<0.05)), pDiff))

%% Figure 1L (histogram of preferred directions)
cols = lines(3);
binSize = 10;
edges = 0:binSize:360;
bins = edges(1:end-1) + binSize/2;
N1 = histcounts(prefDirs(p_DSI < 0.05), edges);
N2 = histcounts(prefDirs(p_DSI >= 0.05 & p_OSI < 0.05), edges);
figure
b = bar(bins, [N2; N1], 'stacked');
b(1).FaceColor = 'w';
b(2).FaceColor = 'k';
hold on
Y = fft(N1+N2);
bins10 = bins + (0:9)'./10.*binSize;
bins10 = bins10(:)';
Y4 = zeros(1,length(bins10));
Y4(1) = Y(1)*10;
Y4([1+4 end-4+1]) = Y([1+4 end-4+1]).*10;
plot(bins10, ifft(Y4), 'r', 'LineWidth', 1)
for ex = 1:size(examplesTun,1)
    idx = strcmp(subjects, examplesTun{ex,1}) & ...
        strcmp(dates, examplesTun{ex,2}) & planes==examplesTun{ex,3} & ...
        ids==examplesTun{ex,4};
    if isnan(prefDirs(idx))
        continue
    end
    plot(prefDirs(idx), 200, 'v', 'Color', cols(ex,:))
end
xlim([0 360])
set(gca, 'XTick', 0:90:360, 'box', 'off')
legend('direction selective','only orientation selective')
xlabel('Preferred direction (deg)')
ylabel('#Boutons')
title(sprintf('n = %d', sum(N1+N2)))