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
exampleEphys = {'SS061', '2016-05-11', 398};

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

%% Prepare for Figs. 4C+D+H
% Collect relevant variables from visual noise data
numSh = 200;
subjects = {};
dates = {};
planes = [];
ids = [];
prefDirs = [];
minima = [];
maxima = [];
means = [];
nullMaxima = [];
nullMeans = [];
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
        pars = cat(3, parsS, parsL);
        prefD = parsS(:,1);
        prefD(isnan(parsS(:,2))) = NaN;
        prefDirs = [prefDirs; prefD];
        nullParsS = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_tuning.nullParametersSmall.npy'));
        nullParsL = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_tuning.nullParametersLarge.npy'));
        nullPars = cat(4, nullParsS, nullParsL);
        curvesS = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_tuning.curvesSmall.npy'));
        curvesL = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_tuning.curvesLarge.npy'));
        curves = cat(3, curvesS, curvesL);
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
        
        mi = NaN(n,2);
        ma = NaN(n,2);
        mn = NaN(n,2);
        nma = NaN(n,2,numSh);
        nmn = NaN(n,2,numSh);
        for iCell = 1:n
            for cond = 1:2
                if ~isnan(pars(iCell,2,cond)) % tuned
                    pd = pars(iCell,1,cond);
                    mn(iCell,cond) = mean(curves(iCell,:,cond));
                    oris = mod(pd + [0 90 180], 360);
                    resp = gratings.orituneWrappedConditions(pars(iCell,:,cond), oris);
                    ma(iCell,cond) = resp(1);
                    if resp(1)-resp(2) < 0 % suppressed by gratings
                        [mi(iCell,cond),ind] = max(resp(2:3));
                    else
                        [mi(iCell,cond),ind] = min(resp(2:3));
                    end
                    ind = ind + 1;
                    
                    resp = NaN(numSh, 3);
                    crv = NaN(numSh, 360);
                    for sh = 1:numSh
                        oris = mod(nullPars(iCell,1,sh,cond) + [0 90 180], 360);
                        resp(sh,:) = gratings.orituneWrappedConditions( ...
                            nullPars(iCell,:,sh,cond), oris);
                        crv(sh,:) = gratings.orituneWrappedConditions(...
                            nullPars(iCell,:,sh,cond), 1:360);
                    end
                    nmi(iCell,cond,:) = resp(:,ind);
                    nma(iCell,cond,:) = resp(:,1);
                    nmn(iCell,cond,:) = mean(crv,2);
                else
                    mn(iCell,cond) = pars(iCell,1,cond);
                    nmn(iCell,cond,:) = nullPars(iCell,1,:,cond);
                end
            end
        end
        minima = [minima; mi];
        maxima = [maxima; ma];
        means = [means; mn];
        nullMaxima = [nullMaxima; nma];
        nullMeans = [nullMeans; nmn];
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

mi = minima(:,1);
ma = maxima(:,1);
j = isnan(mi); % not tuned when pupil small
mi(j) = means(j,1);
ma(j) = means(j,1);

j = ma - mi < 0; % if suppressed by gratings, swap responses to preferred
                 % and non-preferred stimulus because response to pref is
                 % minimum response
m = mi(j);
mi(j) = ma(j);
ma(j) = m;
b_mi = [-flip(2.^(-3:xLim(1))), 0, 2.^(-3:xLim(2))];
b_ma = [-flip(2.^(-3:yLim(1))), 0, 2.^(-3:yLim(2))];

groups = [isGad == -1, isGad == 1];
isTuned = ~any(isnan(maxima),2);
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
    ylabel('#Neurons')
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

%% Figure 4F (example tuning curves, small and large pupil)
cols = 'kr';
for ex = 1:size(examplesTun,1)
    planes = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pRois._ss_2pPlanes.npy'));
    ids = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pRois.ids.npy'));
    directions = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingID.directions.npy'));
    largePupil = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingTrials.largePupil.npy'));
    curvesSmall = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_tuning.curvesSmall.npy'));
    curvesLarge = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_tuning.curvesLarge.npy'));
    curves = cat(3, curvesSmall, curvesLarge);
    amplitudes = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingTrials.amplitudes.npy'));
    
    cellID = examplesTun{ex,3}==planes & examplesTun{ex,4}==ids;
    directions(isnan(directions)) = [];
    
    figure
    hold on
    for cond = 1:2
        amps = amplitudes(:,:,cellID);
        if cond == 1
            ind = largePupil;
        else
            ind = ~largePupil;
        end
        amps(ind) = NaN;
        plot(1:360, curves(cellID,:,cond), 'Color', cols(cond), 'LineWidth',2);
        m = nanmean(amps,2);
        s = nanstd(amps,0,2) ./ sqrt(sum(~ind,2));
        errorbar([directions; 360], m([1:end 1]), s([1:end 1]), 'o', ...
            'Color', cols(cond), 'CapSize', 2, 'MarkerFaceColor', cols(cond))
    end
    plot([0 360], [0 0], 'k:', 'LineWidth', 2)
    set(gca,'XTick',0:90:360)
    title(sprintf('Example %d',ex))
    xlim([0 360])
    xlabel('Direction (in degrees)')
    ylabel('\DeltaF/F')
end

%% Figure 4G (correlations with pupil during gratings)
% Collect data for SC neurons
subjects = {};
dates = {};
planes = [];
ids = [];
rhosNeurons = [];
nullsNeurons = [];
subjDirs = dir(fullfile(folderBase, 'sc neurons 2p', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'sc neurons 2p', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        
        if ~isfile(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsPupil.rhosGratings.npy'))
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
        
        rho = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_corrsPupil.rhosGratings.npy'));
        null = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_corrsPupil.nullRhosGratings.npy'));
        rhosNeurons = [rhosNeurons; rho];
        nullsNeurons = [nullsNeurons; null];
    end
end
pVals = sum(nullsNeurons < rhosNeurons, 2) ./ size(nullsNeurons, 2);
ind = pVals > 0.5;
pVals(ind) = 1 - pVals(ind);
pVals = 2 .* pVals; % two-sided test
pVals(pVals==0) = 1/size(nullsNeurons,2);

% Collect data for retinal boutons
rhosBoutons = [];
subjDirs = dir(fullfile(folderBase, 'boutons', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'boutons', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsPupil.rhosGratings.npy'))
            continue
        end        
        rho = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_corrsPupil.rhosGratings.npy'));
        rhosBoutons = [rhosBoutons; rho];
    end
end

% Make plots
% Histogram
edges = -0.55 : 0.1 : 0.55;
bins = edges(2:end)-0.05;
valid = ~isnan(rhosNeurons);
n1 = histcounts(rhosNeurons(pVals<0.05 & valid),edges)';
n2 = histcounts(rhosNeurons(pVals>=0.05 & valid),edges)';
figure;
b = bar(bins, [n1,n2], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
xlabel(sprintf('Correlation with pupil'))
ylabel('#Neurons')
xlim(edges([1 end]))
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p \geq 0.05')
title(sprintf('n = %d', sum(n1+n2)))

% Cumulative plot
ind = isnan(rhosNeurons);
rhos = rhosNeurons;
rhos(ind) = [];
x = sort(rhos, 'ascend');
y = (1:length(x))' ./ length(x);
x = [-1; x; 1];
y = [0; y; 1];
nulls = nullsNeurons;
nulls(ind,:) = [];
xNull = sort(nulls, 1, 'ascend');
xNull = sort(xNull, 2, 'ascend');
limNull = prctile(xNull, [2.5 97.5], 2);
limNull = [[-1 -1]; limNull; [1 1]];
yNull = (1:length(x))' ./ (length(x));

ind = isnan(rhosBoutons);
rhosB = rhosBoutons;
rhosB(ind) = [];
xB = sort(rhosB, 'ascend');
yB = (1:length(xB))' ./ length(xB);
xB = [-1; xB; 1];
yB = [0; yB; 1];

rhosEx = NaN(size(examplesCorr,1),1);
for ex = 1:size(examplesCorr,1)
    idx = strcmp(subjects, examplesCorr{ex,1}) & ...
        strcmp(dates, examplesCorr{ex,2}) & planes == examplesCorr{ex,3} & ...
        ids == examplesCorr{ex,4};
    rhosEx(ex) = rhosNeurons(idx);
end

cols = 'rb';
figure
plot([0 0], [0 1], 'k')
hold on
h = [0 0 0];
h(3) = plot(xB, yB, 'Color', [.65,.16,.16], 'LineWidth', 2);
h(2) = fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
    'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
h(1) = plot(x, y, 'k', 'LineWidth', 2);
heights = interp1(x, y, rhosEx);
for ex = 1:length(rhosEx)
    plot(rhosEx(ex), heights(ex), 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', cols(ex))
end
xlim([-0.55 0.55])
xlabel('Correlation with pupil')
ylabel('Proportion of neurons')
legend(h, {'sSC neurons','shifted','boutons'}, 'Location', 'NorthWest')
legend('boxoff')
set(gca, 'XTick', [-.5 0 .5], 'box', 'off');

%% Figure 4H (response modulations)
modFun = @(a,b) (b-a)./((abs(a)+abs(b)) ./ 2) .* 100;

% invert sign of responses of suppressed cells
minima(isSuppr==1,:) = -minima(isSuppr==1,:);
maxima(isSuppr==1,:) = -maxima(isSuppr==1,:);
means(isSuppr==1,:) = -means(isSuppr==1,:);
nullMaxima(isSuppr==1,:,:) = -nullMaxima(isSuppr==1,:,:);
nullMeans(isSuppr==1,:,:) = -nullMeans(isSuppr==1,:,:);

mx = maxima;
ind = all(isnan(maxima),2);
mx(ind,:) = means(ind,:);
nmx = nullMaxima;
nmx(ind,:,:) = nullMeans(ind,:,:);
respMod = modFun(mx(:,1), mx(:,2));
nullMod = modFun(squeeze(nmx(:,1,:)), squeeze(nmx(:,2,:)));
confInt = prctile(nullMod, [2.5 97.5], 2);
sgnfcnt = respMod < confInt(:,1) | respMod > confInt(:,2);

% collect data for retinal boutons
maximaB = [];
meansB = [];
isSupprB = [];
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
        parsS = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.parametersSmall.npy'));
        parsL = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.parametersLarge.npy'));
        pars = cat(3, parsS, parsL);
        curvesS = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.curvesSmall.npy'));
        curvesL = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.curvesLarge.npy'));
        curves = cat(3, curvesS, curvesL);
        isSupprB = [isSupprB; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.isSuppressed.npy'))];
        n = size(parsS,1);
        
        ma = NaN(n,2);
        mn = NaN(n,2);
        for iCell = 1:n
            for cond = 1:2
                if ~isnan(pars(iCell,2,cond)) % tuned
                    pd = pars(iCell,1,cond);
                    mn(iCell,cond) = mean(curves(iCell,:,cond));
                    oris = mod(pd + [0 90 180], 360);
                    resp = gratings.orituneWrappedConditions(pars(iCell,:,cond), oris);
                    ma(iCell,cond) = resp(1);
                else
                    mn(iCell,cond) = pars(iCell,1,cond);
                end
            end
        end
        maximaB = [maximaB; ma];
        meansB = [meansB; mn];
    end
end
maximaB(isSupprB==1,:) = -maximaB(isSupprB==1,:);
meansB(isSupprB==1,:) = -meansB(isSupprB==1,:);
mx = maximaB;
ind = all(isnan(maximaB),2);
mx(ind,:) = meansB(ind,:);
respModB = modFun(mx(:,1), mx(:,2));

% Histogram
cols = lines(4);
binSize = 20;
mini = -80;
maxi = 80;
figure
bins = mini:binSize:maxi;
edges = [bins-binSize/2, maxi+binSize/2];
n1 = histcounts(respMod(sgnfcnt), edges);
n2 = histcounts(respMod(~sgnfcnt), edges);
b = bar(bins, [n1',n2'], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
hold on
plot(nanmean(respMod(sgnfcnt & respMod<0)), 700, 'vk', 'MarkerFaceColor', 'k')
plot(nanmean(respMod(sgnfcnt & respMod>0)), 700, 'vk', 'MarkerFaceColor', 'k')
xlim(edges([1 end]))
title(sprintf('n = %d', sum(~isnan(respMod))))
xlabel('Response modulation (%)')
ylabel('#Neurons')
legend(b, 'p < 0.05', 'p \geq 0.05', 'Location', 'NorthWest')
ax = gca;
ax.Box = 'off';
ax.XTick = [mini 0 maxi];

% cumulative distribution
ind = ~isnan(respMod);
x = sort(respMod(ind), 'ascend');
y = (1:sum(ind)) ./ sum(ind);
x = [-200; x; 200];
y = [0 y 1];
pseudo = nullMod;
pseudo(~ind,:) = [];
xNull = sort(pseudo, 1, 'ascend');
xNull = sort(xNull, 2, 'ascend');
limNull = prctile(xNull, [2.5 97.5], 2);
limNull = [[-1 -1]; limNull; [1 1]];
yNull = (1:length(x)) ./ length(x);
ind = ~isnan(respModB);
xB = sort(respModB(ind), 'ascend');
yB = (1:sum(ind)) ./ sum(ind);
xB = [-200; xB; 200];
yB = [0 yB 1];

figure
plot([0 0], [0 1], 'k')
hold on
h = [0 0 0];
h(3) = plot(xB, yB, 'Color', [.65,.16,.16], 'LineWidth', 2);
h(2) = fill([limNull(:,1);flip(limNull(:,2))], [yNull, flip(yNull)], ...
    'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
h(1) = plot(x, y, 'k', 'LineWidth', 2);
[xUni, ind] = unique(x);
yUni = y(ind);
for ex = 1:size(examplesTun,1)
    idx = strcmp(subjects, examplesTun{ex,1}) & ...
        strcmp(dates, examplesTun{ex,2}) & planes == examplesTun{ex,3} & ...
        ids == examplesTun{ex,4};
    md = respMod(idx);
    if md < edges(1)
        md = edges(1);
    end
    height = interp1(xUni, yUni, md);
    if sgnfcnt(idx)
        plot(md, height, 'o', 'MarkerEdgeColor', cols(ex,:), ...
            'MarkerFaceColor', cols(ex,:), 'LineWidth', 2)
    else
        plot(md, height, 'o', 'MarkerEdgeColor', cols(ex,:), 'LineWidth', 2)
    end
end
xlim(edges([1 end]))
xlabel('Response modulation (%)')
ylabel('Proportion of neurons')
title(sprintf('n = %d', sum(~isnan(respMod))))
legend(h, {'sSC neurons','shifted','boutons'}, 'Location', 'NorthWest')
legend('boxoff')
set(gca, 'XTick', [mini 0 maxi], 'box', 'off');

%% Figure 4I (ephys: example tuning curves, small/large pupil, with/without laser)
cols = 'kr';

ids = readmatrix(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\clusters.uuids.csv'));
directions = readNPY(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\_ss_gratingID.directions.npy'));
laserOn = readNPY(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\_ss_gratingID.laserOn.npy'));
largePupil = readNPY(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\_ss_gratingTrials.largePupil.npy'));
curvesSmallLaserOff = readNPY(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\_ss_tuning.curvesSmallLaserOff.npy'));
curvesLargeLaserOff = readNPY(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\_ss_tuning.curvesLargeLaserOff.npy'));
curvesSmallLaserOn = readNPY(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\_ss_tuning.curvesSmallLaserOn.npy'));
curvesLargeLaserOn = readNPY(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\_ss_tuning.curvesLargeLaserOn.npy'));
curves = cat(3, curvesSmallLaserOff, curvesLargeLaserOff, ...
    curvesSmallLaserOn, curvesLargeLaserOn);
amplitudes = readNPY(fullfile(folderBase, 'sc neurons ephys', exampleEphys{1}, ...
    exampleEphys{2}, '001\_ss_gratingTrials.amplitudes.npy'));
    
cellID = exampleEphys{3}==ids;
laserConds = [false false; true true];
pupilConds = [false true; false true];
ttls = {'Control', 'V1 inactivated'};
for ls = 1:2
    figure
    hold on
    for cond = 1:2
        amps = amplitudes(:,:,cellID);
        ind = largePupil==pupilConds(ls,cond) & laserOn==laserConds(ls,cond);
        amps(~ind) = NaN;
        plot(1:360, curves(cellID,:,(ls-1)*2+cond), 'Color', cols(cond), 'LineWidth',2);
        ind2 = ~all(isnan(amps),2) & ~isnan(directions);
        m = nanmean(amps(ind2,:),2);
        s = nanstd(amps(ind2,:),0,2) ./ sqrt(sum(ind(ind2,:),2));
        errorbar([directions(ind2); 360], m([1:end 1]), s([1:end 1]), 'o', ...
            'Color', cols(cond), 'CapSize', 2, 'MarkerFaceColor', cols(cond))
        plot([0 360], [1 1].*nanmean(amps(~all(isnan(amps),2) & ...
            isnan(directions),:),2), ['--' cols(cond)], 'LineWidth', 2)
    end
    set(gca,'XTick',0:90:360)
    title(ttls{ls})
    xlim([0 360])
    ylim([-3 35])
    xlabel('Direction (in degrees)')
    ylabel('Firing rate (spikes/s)')
end

%% Prepare for Figs. J+K
numSh = 200;
subjects = {};
dates = {};
ids = [];
maxima = [];
means = [];
nullBehMaxima = [];
nullBehMeans = [];
nullLaserMaxima = [];
nullLaserMeans = [];
isSuppr = [];

subjDirs = dir(fullfile(folderBase, 'sc neurons ephys'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    if ~subjDirs(subj).isdir || any(strcmp(name, {'.','..'}))
        continue
    end
    dateDirs = dir(fullfile(folderBase, 'sc neurons ephys', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        c = readmatrix(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\clusters.uuids.csv'));
        n = length(c);
        ids = [ids; c];
        subjects = [subjects; repmat({name}, n, 1)];
        dates = [dates; repmat({date}, n, 1)];
        
        parsSOff = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.parametersSmallLaserOff.npy'));
        parsLOff = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.parametersLargeLaserOff.npy'));
        parsSOn = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.parametersSmallLaserOn.npy'));
        parsLOn = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.parametersLargeLaserOn.npy'));
        pars = cat(3, parsSOff, parsLOff, parsSOn, parsLOn);
        cSOff = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.curvesSmallLaserOff.npy'));
        cLOff = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.curvesLargeLaserOff.npy'));
        cSOn = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.curvesSmallLaserOn.npy'));
        cLOn = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.curvesLargeLaserOn.npy'));
        curves = cat(3, cSOff, cLOff, cSOn, cLOn);
        nbSOff = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.nullBehaviourParametersSmallLaserOff.npy'));
        nbLOff = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.nullBehaviourParametersLargeLaserOff.npy'));
        nbSOn = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.nullBehaviourParametersSmallLaserOn.npy'));
        nbLOn = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.nullBehaviourParametersLargeLaserOn.npy'));
        nbPars = cat(4, nbSOff, nbLOff, nbSOn, nbLOn);
        nlSOff = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.nullLaserParametersSmallLaserOff.npy'));
        nlLOff = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.nullLaserParametersLargeLaserOff.npy'));
        nlSOn = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.nullLaserParametersSmallLaserOn.npy'));
        nlLOn = readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.nullLaserParametersLargeLaserOn.npy'));
        nlPars = cat(4, nlSOff, nlLOff, nlSOn, nlLOn);
        isSuppr = [isSuppr; readNPY(fullfile(folderBase, 'sc neurons ephys', name, ...
            date, '001\_ss_tuning.isSuppressed.npy'))];
        
        ma = NaN(n,4);
        mn = NaN(n,4);
        nbma = NaN(n,4,numSh);
        nbmn = NaN(n,4,numSh);
        nlma = NaN(n,4,numSh);
        nlmn = NaN(n,4,numSh);
        for iCell = 1:n
            for cond = 1:4
                if ~isnan(pars(iCell,2,cond)) % tuned
                    pd = pars(iCell,1,cond);
                    mn(iCell,cond) = mean(curves(iCell,:,cond));
                    oris = mod(pd + [0 90 180], 360);
                    resp = gratings.orituneWrappedConditions(pars(iCell,:,cond), oris);
                    ma(iCell,cond) = resp(1);
                    
                    respB = NaN(numSh, 3);
                    crvB = NaN(numSh, 360);
                    respL = NaN(numSh, 3);
                    crvL = NaN(numSh, 360);
                    for sh = 1:numSh
                        oris = mod(nbPars(iCell,1,sh,cond) + [0 90 180], 360);
                        respB(sh,:) = gratings.orituneWrappedConditions( ...
                            nbPars(iCell,:,sh,cond), oris);
                        crvB(sh,:) = gratings.orituneWrappedConditions(...
                            nbPars(iCell,:,sh,cond), 1:360);
                        oris = mod(nlPars(iCell,1,sh,cond) + [0 90 180], 360);
                        respL(sh,:) = gratings.orituneWrappedConditions( ...
                            nlPars(iCell,:,sh,cond), oris);
                        crvL(sh,:) = gratings.orituneWrappedConditions(...
                            nlPars(iCell,:,sh,cond), 1:360);
                    end
                    nbma(iCell,cond,:) = respB(:,1);
                    nbmn(iCell,cond,:) = mean(crvB,2);
                    nlma(iCell,cond,:) = respL(:,1);
                    nlmn(iCell,cond,:) = mean(crvL,2);
                else
                    mn(iCell,cond) = pars(iCell,1,cond);
                    nbmn(iCell,cond,:) = nbPars(iCell,1,:,cond);
                    nlmn(iCell,cond,:) = nlPars(iCell,1,:,cond);
                end
            end
        end
        maxima = [maxima; ma];
        means = [means; mn];
        nullBehMaxima = [nullBehMaxima; nbma];
        nullBehMeans = [nullBehMeans; nbmn];
        nullLaserMaxima = [nullLaserMaxima; nlma];
        nullLaserMeans = [nullLaserMeans; nlmn];
    end
end

% for untuned neurons, set maxima to mean responses
isTuned = ~isnan(maxima(:,1));
maxima(~isTuned,:) = means(~isTuned,:);
nullBehMaxima(~isTuned,:,:) = nullBehMeans(~isTuned,:,:);
nullLaserMaxima(~isTuned,:,:) = nullLaserMeans(~isTuned,:,:);

modFun = @(a,b) (b-a)./((abs(a)+abs(b)) ./ 2) .* 100;

% response modulation during control and V1 inactivation
diffOff = modFun(maxima(:,1), maxima(:,2));
diffOn = modFun(maxima(:,3), maxima(:,4));
validUnits = ~isnan(diffOff) & ~isnan(diffOn);
% response modulation during control and V1 inactivation for surrogate data
% where label of pupil size has been randomised
pseudoBehOff = modFun(squeeze(nullBehMaxima(:,1,:)), squeeze(nullBehMaxima(:,2,:)));
pseudoBehOn = modFun(squeeze(nullBehMaxima(:,3,:)), squeeze(nullBehMaxima(:,4,:)));
% response modulation during control and V1 inactivation for surrogate data
% where label of laser on/off has been randomised
pseudoLaserOff = modFun(squeeze(nullLaserMaxima(:,1,:)), squeeze(nullLaserMaxima(:,2,:)));
pseudoLaserOn = modFun(squeeze(nullLaserMaxima(:,3,:)), squeeze(nullLaserMaxima(:,4,:)));

% for suppressed units, change sign of response modulations so that
% a more negative response with arousal results in a positive DI
diffOff(isSuppr==1) = -diffOff(isSuppr==1);
diffOn(isSuppr==1) = -diffOn(isSuppr==1);
pseudoBehOff(isSuppr==1,:) = -pseudoBehOff(isSuppr==1,:);
pseudoBehOn(isSuppr==1,:) = -pseudoBehOn(isSuppr==1,:);
pseudoLaserOff(isSuppr==1,:) = -pseudoLaserOff(isSuppr==1,:);
pseudoLaserOn(isSuppr==1,:) = -pseudoLaserOn(isSuppr==1,:);

% significance of response modulation during control condition (using
% surrogate data where pupil size label was randomised)
confIntBehOff = prctile(pseudoBehOff, [2.5 97.5], 2);
signBehOff = diffOff < confIntBehOff(:,1) | diffOff > confIntBehOff(:,2);
% significance of response modulation during V1 inactivation (using
% surrogate data where pupil size label was randomised)
confIntBehOn = prctile(pseudoBehOn, [2.5 97.5], 2);
signBehOn = diffOn < confIntBehOn(:,1) | diffOn > confIntBehOn(:,2);

% difference between response modulation during control vs V1 inactivation
% + significance of the difference (using surrogate data where laser label
% has been randomised)
diffs = diffOn - diffOff;
pseudoDiffsLaser = pseudoLaserOn - pseudoLaserOff;
confIntLaser = prctile(pseudoDiffsLaser, [2.5 97.5], 2);
signLaser = diffs < confIntLaser(:,1) | diffs > confIntLaser(:,2);

%% Figure 4J (ephys: response modulation histograms)
binSizes = 20;
bins = -200 : binSizes : 200;
edges = [bins bins(end)+binSizes] - binSizes/2;

% histogram for laser off condition
figure
n1 = histcounts(diffOff(signBehOff & validUnits), edges);
n2 = histcounts(diffOff(~signBehOff & validUnits), edges);
b = bar(bins, [n1',n2'], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
xlim([-210 210])
ylim([0 70])
title(sprintf('Control (n = %d)', sum(validUnits)))
xlabel('Response modulation (%)')
ylabel('#Neurons')
set(gca, 'box', 'off', 'XTick', [-200 0 200])

% histogram for laser on condition
figure
n1 = histcounts(diffOn(signBehOn & validUnits), edges);
n2 = histcounts(diffOn(~signBehOn & validUnits), edges);
b = bar(bins, [n1',n2'], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
xlim([-210 210])
ylim([0 70])
title(sprintf('V1 inactivated (n = %d)', sum(validUnits)))
xlabel('Response modulation (%)')
ylabel('#Neurons')
set(gca, 'box', 'off', 'XTick', [-200 0 200])
    
%% Figure 4K (ephys: scatter of response modulation during control vs V1 inactivation)

n = [0 0];
h = [0 0];
figure('Position', [540 560 700 420])
hold on
% (1) laser has no significant impact on response modulation (large dots) (small dots)
ind = ~signLaser & validUnits;
n(1) = sum(ind);
h(1) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 7, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1] .* .5);
% (2) laser has significant impact on response modulation (large dots)
ind = signLaser & validUnits;
n(2) = sum(ind);
h(2) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 7, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
exID = strcmp(subjects, exampleEphys{1}) & strcmp(dates, exampleEphys{2}) & ...
    ids==exampleEphys{3};
plot(diffOff(exID), diffOn(exID), 'o', ...
    'MarkerSize', 5, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r')
    
plot([-200 200], [-200 200], 'k:', 'LineWidth', 1)
legend(h, ['p_{laser} \geq 0.05 (n=' num2str(n(1)) ')'], ...
    ['p_{laser} < 0.05 (n=' num2str(n(2)) ')'], ...
    'Location', 'northeastoutside')
axis([-1 1 -1 1].*200)
axis square
xlabel(sprintf('Response modulation (%%)\ncontrol'))
ylabel(sprintf('Response modulation (%%)\nV1 inactivation'))
title(sprintf('n = %d', sum(n)))
set(gca, 'box', 'off', 'XTick', [-200 0 200], 'YTick', [-200 0 200])