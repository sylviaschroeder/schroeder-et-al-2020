%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish';
folderTools = 'C:\STORAGE\workspaces';
folderThisRepo = 'C:\dev\workspace\schroeder-et-al-2020';

%% Parameters
% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

% Preprocess neural data
sigma = 1; % in sec

% to group into ON, OFF and ON+OFF
onThr = 0.5; % ON field is at least three times stronger than OFF field
offThr = -0.5; % OFF field is at least three times stronger than ON field

%% Examples
examplesTun = {'SS038', '2015-02-17', 1, 227; ...
               'SS041', '2015-04-23', 2, 38; ...
               'SS041', '2015-04-23', 2, 151; ...
               'SS038', '2015-02-17', 1, 203};
examplesCorr = {'SS047', '2015-11-23', 1, 153; ...
                'SS047', '2015-11-23', 1, 171};

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Figure 4B (histogram of ON/OFF indices)
% Collect relevant variables from visual noise data
subjects = {};
dates = {};
planes = [];
ids = [];
evStim = [];
lambdasStim = [];
pValues = [];
OnOffValues = [];

subjDirs = dir(fullfile(folderBase, 'sc neurons 2p', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'sc neurons 2p', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_rf.maps.npy'))
            continue
        end
        
        planes = [planes; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_2pRois._ss_2pPlanes.npy'))];
        ids = [ids; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_2pRois.ids.npy'))];
        rfs_exp = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_rf.maps.npy'));
        evStim = [evStim; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_rf.explVarsStim.npy'))];
        lambdasStim = [lambdasStim; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_rf.lambdasStim.npy'))];
        pValues = [pValues; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
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
set(gca, 'box', 'off', 'XTick', [-1 offThr 0 onThr 1])
xlim([-1 1])
ylim([0 mx])
xlabel('ON/OFF index')
ylabel('#sc neurons 2p')
title(sprintf('n = %d', sum(validRF)))

%% Prepare for Figs. 4C+D
% Collect relevant variables from visual noise data
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
isGad = [];

subjDirs = dir(fullfile(folderBase, 'sc neurons 2p', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'sc neurons 2p', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_tuning.parametersSmall.npy'))
            continue
        end
        
        p = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_2pRois._ss_2pPlanes.npy'));
        n = length(p);
        planes = [planes; p];
        ids = [ids; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_2pRois.ids.npy'))];
        subjects = [subjects; repmat({name}, n, 1)];
        dates = [dates; repmat({date}, n, 1)];
        parsS = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_tuning.parametersSmall.npy'));
        parsL = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_tuning.parametersLarge.npy'));
        isT = [~isnan(parsS(:,2)), ~isnan(parsL(:,2))];
        isTuned = [isTuned; all(isT,2)];
        prefD = parsS(:,1);
        prefD(~isT(:,1)) = NaN;
        prefDirs = [prefDirs; prefD];
        curves = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_tuning.curvesSmall.npy'));
        isSuppr = [isSuppr; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_tuning.isSuppressed.npy'))];
        isGad = [isGad; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_2pRois.isGad.npy'))];
        amp = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_gratingTrials.amplitudes.npy'));
        amplitudes = [amplitudes; permute(mat2cell(amp, ...
            size(amp,1), size(amp,2), ones(1,n)), [3 1 2])];
        largePupil = [largePupil; repmat({readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_gratingTrials.largePupil.npy'))}, n, 1)];
        
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

%% Figure 4C (histograms of max responses for exc. and inh. neurons)
cols = lines(4);
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

groups = [isGad == -1, isGad == 1];

for g = 1:2
    figure
    N1 = histcounts(ma(isTuned==1 & groups(:,g)), b_ma);
    N2 = histcounts(ma(isTuned==0 & groups(:,g)), b_ma);
    b = bar(1.5:length(b_ma), [N1',N2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    hold on
    plot([8 8], [0 max(N1+N2)*1.1], 'k')
    for ex = 1:size(examplesTun,1)
        idx = strcmp(subjects, examplesTun{ex,1}) & ...
            strcmp(dates, examplesTun{ex,2}) & planes==examplesTun{ex,3} & ...
            ids==examplesTun{ex,4};
        if ~groups(idx,g)
            continue
        end
        exMax = find(b_ma<ma(idx),1,'last');
        if ma(idx)>0
            exMax = exMax + rem(log2(ma(idx)),1);
        else
            exMax = exMax - rem(log2(-ma(idx)),1);
        end
        plot(exMax,max(N1+N2)*1.04,'v', 'MarkerFaceColor', cols(ex,:), ...
            'MarkerEdgeColor', 'none')
    end
    ylim([0 max(N1+N2)*1.1])
    xlabel('Maximum')
    ylabel('#sc neurons 2p')
    title(sprintf('n = %d', sum(N1 + N2)))
    set(gca,'XTick',[1 4 8 12 15 17],'XTickLabel',b_ma([1 4 8 12 15 17]),'box','off')
    if g == 1
        legend(b,{'tuned','not tuned'},'Location','NorthWest')
    end
end

%% Figure 4D (histogram of preferred directions)
cols = lines(4);
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
    plot(prefDirs(idx), 80, 'v', 'MarkerFaceColor', cols(ex,:), ...
        'MarkerEdgeColor', 'none')
end
xlim([0 360])
set(gca, 'XTick', 0:90:360, 'box', 'off')
legend('direction selective','only orientation selective', 'Location', 'NorthWest')
xlabel('Preferred direction (deg)')
ylabel('#Neurons')
title(sprintf('n = %d', sum(N1+N2)))

%% Figure 4E (neuron traces during gratings)
planes = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pRois._ss_2pPlanes.npy'));
ids = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pRois.ids.npy'));
traces = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pCalcium.dff.npy'));
time = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pCalcium.timestamps.npy'));
delays = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pPlanes.delay.npy'));
running = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_running.speed.npy'));
runningTime = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_running.timestamps.npy'));
pupil = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\eye.diameter.npy'));
pupilTime = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\eye.timestamps.npy'));
gratingTimes = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_grating.intervals.npy'));
recTimeGratings = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_recordings.gratings_intervals.npy'));

ind = time >= recTimeGratings(1) & time <= recTimeGratings(2);
traces = traces(ind,:);
time = time(ind);

frameDur = median(diff(time));
stdSamples = round(sigma / frameDur);
convWindow = normpdf(-4*stdSamples:4*stdSamples, 0, stdSamples);

ind = isnan(running);
indInterp = histcounts(runningTime(ind), time) > 0;
runningNew = interp1(runningTime(~ind), running(~ind), time, 'pchip');
runningNew = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
    mean(runningNew(1:round(length(convWindow)/2))); ...
    runningNew; ones(floor((length(convWindow)-1)/2),1) .* ...
    mean(runningNew(end-round(length(convWindow)/2):end))], ...
    convWindow, 'valid');
runningNew(indInterp) = NaN;
ind = isnan(pupil);
indInterp = histcounts(pupilTime(ind), time) > 0;
pupilNew = interp1(pupilTime(~ind), pupil(~ind), time, 'pchip');
pupilNew = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
    mean(pupilNew(1:round(length(convWindow)/2))); ...
    pupilNew; ones(floor((length(convWindow)-1)/2),1) .* ...
    mean(pupilNew(end-round(length(convWindow)/2):end))], ...
    convWindow, 'valid');
pupilNew(indInterp) = NaN;

filteredTraces = NaN(size(traces));
for n = 1:size(traces,2)
    tr = traces(:,n);
    ind = isnan(tr);
    tr(ind) = interp1(time(~ind), tr(~ind), time(ind), 'pchip');
    filteredTraces(:,n) = ...
        conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
        mean(tr(1:round(length(convWindow)/2))); ...
        tr; ...
        ones(floor((length(convWindow)-1)/2),1) .* ...
        mean(tr(end-round(length(convWindow)/2)))], ...
        convWindow, 'valid');
end
centerIndex = floor(size(traces,1)/2);
ind = true(size(traces,1),1);
ind(centerIndex+1:end) = false;
ind(isnan(runningNew)) = false;
rhos = corr(runningNew(ind), filteredTraces(ind,:))';

[r_p_sorted,order] = sort(rhos,'descend');
ids_sorted = ids(order,:);
zscored = (filteredTraces - nanmean(filteredTraces,1)) ./ ...
    nanstd(filteredTraces,0,1) ./ size(filteredTraces,1).^0.5;
zscored = zscored(:,order);
if ~all(isnan(pupilNew))
    ind = ~isnan(pupilNew);
    pupilNew = interp1(time(ind), pupilNew(ind), time, 'pchip');
end

cols = 'rb';
cm = flip(gray);
ax = [0 0 0];
figure('Position',[15 45 1530 940])
subplot(5,1,1)
hold on
h = [0 0 0];
h(1) = plot(time, pupilNew-20, 'Color', [0 0.7 0.5]);
h(2) = plot(time, runningNew./2, 'Color', [0 0 .7]);
numStim = size(gratingTimes,1);
stimT = [repmat(gratingTimes(:,1)',2,1);NaN(1,numStim)];
s = [zeros(1,numStim); ones(1,numStim).*5; NaN(1,numStim)];
h_ = plot(stimT,s-12,'k');
h(3) = h_(1);
leg = legend(h, 'Running','Pupil','Stimulus');
title('Gratings')
ylim([-20 45])
set(gca,'box','off')
ax(1) = gca;
leg.Position = [0.92,0.83,0.05,0.09];
subplot(5,1,2:4)
imagesc(time([1 end]),[1 size(zscored,2)], zscored', ...
    prctile(zscored(:),[5 95]))
colormap(cm)
set(gca,'box','off')
ax(2) = gca;
subplot(5,1,5)
hold on
for ex = 1:size(examplesCorr,1)
    idx = planes == examplesCorr{ex,3} & ids == examplesCorr{ex,4};
    plot(time, smooth(traces(:,idx),5)-10*ex, cols(ex))
end
ax(3) = gca;
linkaxes(ax,'x')
xlim(time([centerIndex+1 end]))
xlabel('Time (s)')
set(gca,'box','off')