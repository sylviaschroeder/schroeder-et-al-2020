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

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(genpath(fullfile(folderTools, 'CircStat2012a')))
addpath(fullfile(folderThisRepo))

%% Figure S2A (probability of running given pupil size)
% collect running and pupil data during gratings and gray screens from all
% datasets
running = {};
timeRun = {};
pupil = {};
timePupil = {};
timeCa = {};
intervalGratings = [];
intervalGray = [];

subjDirs = dir(fullfile(folderBase, 'sc neurons 2p', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'sc neurons 2p', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        folder = fullfile(folderBase, 'sc neurons 2p', name, date, '001');
        if ~isfile(fullfile(folder, 'eye.diameter.npy')) || ...
                ~isfile(fullfile(folder, '_ss_running.speed.npy'))
            continue
        end
        running{end+1,1} = readNPY(fullfile(folder, '_ss_running.speed.npy'));
        timeRun{end+1,1} = readNPY(fullfile(folder, '_ss_running.timestamps.npy'));
        pupil{end+1,1} = readNPY(fullfile(folder, 'eye.diameter.npy'));
        timePupil{end+1,1} = readNPY(fullfile(folder, 'eye.timestamps.npy'));
        timeCa{end+1,1} = readNPY(fullfile(folder, '_ss_2pCalcium.timestamps.npy'));
        if isfile(fullfile(folder, '_ss_recordings.gratings_intervals.npy'))
            intervalGratings(end+1,:) = readNPY(fullfile(folder, ...
                '_ss_recordings.gratings_intervals.npy'));
        else
            intervalGratings(end+1,:) = NaN(1,2);
        end
        if isfile(fullfile(folder, '_ss_recordings.grayScreen_intervals.npy'))
            intervalGray(end+1,:) = readNPY(fullfile(folder, ...
                '_ss_recordings.grayScreen_intervals.npy'));
        else
            intervalGray(end+1,:) = NaN(1,2);
        end
    end
end

smoothStd = 1; % in sec (to smooth running and pupil traces)
binSize = 0.2; % in sec (to resample running and pupil traces consistently)
intervals = cat(3, intervalGratings, intervalGray);
pupilAll = [];
runningAll = [];
for k = 1:length(running)
    % Gaussian window for convolution of running and pupil traces
    stdSamples = round(smoothStd / median(diff(timeCa{k})));
    convWindow = normpdf(-4*stdSamples:4*stdSamples, 0, stdSamples);
    for j = 1:2 % gratings, and gray screens
        if isnan(intervals(k,1,j))
            continue
        end
        t = timeCa{k}; % time samples of calcium traces
        t = t(find(timeCa{k} <= intervals(k,1,j), 1, 'last') : ...
            find(timeCa{k} >= intervals(k,2,j), 1, 'first')); % only during gratings/gray screens
        tNew = (t(1) : binSize : t(end))'; % new resampled time
        
        r = running{k};
        indsRun = [find(timeRun{k} <= intervals(k,1,j), 1, 'last'), ...
            find(timeRun{k} >= intervals(k,2,j), 1, 'first')];
        r = r(indsRun(1):indsRun(2)); % running only during gratings/gray screens
        rt = timeRun{k}(indsRun(1):indsRun(2));
        % downsample running trace to time of calcium traces
        ind = isnan(r);
        indInterp = histcounts(rt(ind), t) > 0;
        r = interp1(rt(~ind), r(~ind), t, 'pchip');
        % convolve running trace with Gaussian window
        r = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
            mean(r(1:round(length(convWindow)/2))); ...
            r; ...
            ones(floor((length(convWindow)-1)/2),1) .* ...
            mean(r(end-round(length(convWindow)/2):end))], ...
            convWindow, 'valid');
        % downsample running trace to new time scale
        r = interp1(t(~indInterp), r(~indInterp), tNew, 'pchip');
        % normalize trace so all values are between 0 and 1
        r = (r - min(r)) ./ (max(r) - min(r));
        runningAll = [runningAll; r];
        
        p = pupil{k}; % do the same with pupil trace
        indsPup = [find(timePupil{k} <= intervals(k,1,j), 1, 'last'), ...
            find(timePupil{k} >= intervals(k,2,j), 1, 'first')];
        p = p(indsPup(1):indsPup(2));
        pt = timePupil{k}(indsPup(1):indsPup(2));
        ind = isnan(p);
        indInterp = histcounts(pt(ind), t) > 0;
        p = interp1(pt(~ind), p(~ind), t, 'pchip');
        p = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
            mean(p(1:round(length(convWindow)/2))); ...
            p; ...
            ones(floor((length(convWindow)-1)/2),1) .* ...
            mean(p(end-round(length(convWindow)/2):end))], ...
            convWindow, 'valid');
        p = interp1(t(~indInterp), p(~indInterp), tNew, 'pchip');
        p = (p - min(p)) ./ (max(p) - min(p));
        pupilAll = [pupilAll; p];
    end
end

% make conditional probability plot
pBins = prctile(pupilAll,0:10:100); % percentiles of pupil data define bin edges of plot (x-axis)
pBins(end) = pBins(end) + 1;
rBins = prctile(runningAll,0:10:100); % percentiles of pupil data define bin edges of plot (y-axis)
prob = NaN(length(rBins)-1, length(pBins)-1);
for b = 1:length(pBins)-1
    ind = pupilAll>=pBins(b) & pupilAll<pBins(b+1);
    n = histcounts(runningAll(ind), rBins);
    prob(:,b) = n ./ sum(ind);
end
figure
imagesc([5 95], [5 95], prob)
h = colorbar;
h.Label.String = 'P(running|pupil)';
colormap(flip(gray))
ax = gca;
ax.YDir = 'normal';
ax.Box = 'off';
axis square
xlabel('Pupil size (percentile)')
ylabel('Running speed (percentile)')

%% Prepare for rest of plots
numSh = 200;
subjects = {};
dates = {};
planes = [];
ids = [];

minima = [];
maxima = [];
means = [];
nullMinima = [];
nullMaxima = [];
nullMeans = [];
isSuppr = [];
amplitudes = {};
largePupil = {};
directions = {};

rhosRunDark = [];
nullsRunDark = [];
rhosRunGratings = [];
nullsRunGratings = [];
rhosPupilGratings = [];
nullsPupilGratings = [];

evStim = [];
lambdasStim = [];
pValues = [];
OnOffValues = [];

subjDirs = dir(fullfile(folderBase, 'boutons', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    fprintf('Subject %s (%d of %d)\n', name, subj, length(subjDirs))
    dateDirs = dir(fullfile(folderBase, 'boutons', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        
        p = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_2pRois._ss_2pPlanes.npy'));
        n = length(p);
        planes = [planes; p];
        ids = [ids; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_2pRois.ids.npy'))];
        subjects = [subjects; repmat({name}, n, 1)];
        dates = [dates; repmat({date}, n, 1)];
        
        if isfile(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsRunning.rhosDark.npy'))
            rho = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsRunning.rhosDark.npy'));
            null = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsRunning.nullRhosDark.npy'));
        else
            rho = NaN(n,1);
            null = NaN(n,500);
        end
        rhosRunDark = [rhosRunDark; rho];
        nullsRunDark = [nullsRunDark; null];
        if isfile(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsRunning.rhosGratings.npy'))
            rho = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsRunning.rhosGratings.npy'));
            null = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsRunning.nullRhosGratings.npy'));
        else
            rho = NaN(n,1);
            null = NaN(n,500);
        end
        rhosRunGratings = [rhosRunGratings; rho];
        nullsRunGratings = [nullsRunGratings; null];
        if isfile(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsPupil.rhosGratings.npy'))
            rho = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsPupil.rhosGratings.npy'));
            null = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_corrsPupil.nullRhosGratings.npy'));
        else
            rho = NaN(n,1);
            null = NaN(n,500);
        end
        rhosPupilGratings = [rhosPupilGratings; rho];
        nullsPupilGratings = [nullsPupilGratings; null];
        
        if isfile(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.parametersSmall.npy'))
            parsS = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.parametersSmall.npy'));
            parsL = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.parametersLarge.npy'));
            pars = cat(3, parsS, parsL);
            nullParsS = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.nullParametersSmall.npy'));
            nullParsL = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.nullParametersLarge.npy'));
            nullPars = cat(4, nullParsS, nullParsL);
            curvesS = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.curvesSmall.npy'));
            curvesL = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.curvesLarge.npy'));
            curves = cat(3, curvesS, curvesL);
            isSuppr = [isSuppr; readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_tuning.isSuppressed.npy'))];
            amp = readNPY(fullfile(folderBase, ...
                'boutons', name, date, '001\_ss_gratingTrials.amplitudes.npy'));
            amplitudes = [amplitudes; permute(mat2cell(amp, ...
                size(amp,1), size(amp,2), ones(1,n)), [3 1 2])];
            largePupil = [largePupil; repmat({readNPY(fullfile(folderBase, ...
                'boutons', name, date, '001\_ss_gratingTrials.largePupil.npy'))}, n, 1)];
            directions = [directions; repmat({readNPY(fullfile(folderBase, ...
                'boutons', name, date, '001\_ss_gratingID.directions.npy'))}, n, 1)];
            
            mi = NaN(n,2);
            ma = NaN(n,2);
            mn = NaN(n,2);
            nmi = NaN(n,2,numSh);
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
            nullMinima = [nullMinima; nmi];
            nullMaxima = [nullMaxima; nma];
            nullMeans = [nullMeans; nmn];
        else
            isSuppr = [isSuppr; NaN(n,1)];
            amplitudes = [amplitudes; cell(n,1)];
            largePupil = [largePupil; cell(n,1)];
            directions = [directions; cell(n,1)];
            minima = [minima; NaN(n,2)];
            maxima = [maxima; NaN(n,2)];
            means = [means; NaN(n,2)];
            nullMinima = [nullMinima; NaN(n,2,numSh)];
            nullMaxima = [nullMaxima; NaN(n,2,numSh)];
            nullMeans = [nullMeans; NaN(n,2,numSh)];
        end
        
        if isfile(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_rf.maps.npy'))
            rfs = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_rf.maps.npy'));
            ev = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_rf.explVarsStim.npy'));
            lam = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_rf.lambdasStim.npy'));
            pv = readNPY(fullfile(folderBase, 'boutons', name, ...
                date, '001\_ss_rf.pValues.npy'));
            
            oov = NaN(n,2);
            for iCell = 1:n
                rf = squeeze(rfs(iCell,:,:,:,:));
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
        else
            ev = NaN(n,1);
            lam = NaN(n,1);
            pv = NaN(n,1);
            oov = NaN(n,2);
        end
        evStim = [evStim; ev];
        lambdasStim = [lambdasStim; lam];
        pValues = [pValues; pv];
        OnOffValues = [OnOffValues; oov];
    end
end

% invert sign of responses of suppressed cells
minima(isSuppr==1,:) = -minima(isSuppr==1,:);
maxima(isSuppr==1,:) = -maxima(isSuppr==1,:);
means(isSuppr==1,:) = -means(isSuppr==1,:);
nullMinima(isSuppr==1,:,:) = -nullMinima(isSuppr==1,:,:);
nullMaxima(isSuppr==1,:,:) = -nullMaxima(isSuppr==1,:,:);
nullMeans(isSuppr==1,:,:) = -nullMeans(isSuppr==1,:,:);

% find unique datasets
[~,~,sbj] = unique(subjects);
[~,~,dt] = unique(dates);
[~,indDatasets,dataset] = unique([sbj, dt], 'rows');

% Determine ON/OFF indices, and which boutons have a valid RF
[~,type] = max(abs(OnOffValues),[],2);
ind = sub2ind(size(OnOffValues), (1:size(OnOffValues,1))', type);
signs = sign(OnOffValues(ind));
OnOffRatios = OnOffValues;
OnOffRatios(signs<0,:) = -OnOffRatios(signs<0,:);
OnOffRatios(OnOffRatios<0) = 0;
OnOffRatios = (OnOffRatios(:,1)-OnOffRatios(:,2))./sum(OnOffRatios,2);
validRF = pValues < minPVal & evStim > minExplainedVarianceStim & ...
    lambdasStim < maxLambda;

% Determine DSIsSmall and OSIsSmall
shuffles = 1000;
dirct = 0:30:330;
dirVectors = exp(dirct./180.*pi .* 1i);
oriVectors = exp(dirct./180.*2.*pi .* 1i);
ampSmall = amplitudes; % {nROIs x 1}, each entry: [nStim x nReps]
ampLarge = amplitudes;
ampShuffleSmall = cell(size(amplitudes));  % {nROIs x 1}, each entry: [nReps x nStim x nShuffles]
ampShuffleLarge = cell(size(amplitudes));
sz = cell2mat(cellfun(@size, amplitudes, 'UniformOutput', false));
numReps = setdiff(unique(sz(:,2)),0);
numStim = setdiff(unique(sz(:,1)),0);
permutations = cell(1, length(numReps));
for m = 1:length(numReps)
    permutations{m} = NaN(numReps(m)*numStim,shuffles);
    for sh = 1:shuffles
        permutations{m}(:,sh) = randperm(numReps(m)*numStim);
    end
end
for n = 1:size(amplitudes,1)
    m = find(numReps == size(amplitudes{n},2));
    if all(isnan(amplitudes{n}(:))) || isempty(amplitudes{n})
        ampSmall{n} = NaN(numStim, numReps(m));
        ampLarge{n} = NaN(numStim, numReps(m));
        ampShuffleSmall{n} = NaN(numStim, numReps(m), shuffles);
        ampShuffleLarge{n} = NaN(numStim, numReps(m), shuffles);
        continue
    end
    ampSmall{n}(largePupil{n}) = NaN;
    ampLarge{n}(~largePupil{n}) = NaN;
    a = reshape(ampSmall{n}, [], 1);
    a = a(permutations{m});
    ampShuffleSmall{n} = reshape(a, numStim, numReps(m), shuffles);
    a = reshape(ampLarge{n}, [], 1);
    a = a(permutations{m});
    ampShuffleLarge{n} = reshape(a, numStim, numReps(m), shuffles);
end
meanAmpSmall = cellfun(@nanmean, ampSmall, repmat({2},size(amplitudes,1),1), ...
    'UniformOutput', false);
meanAmpSmall = cell2mat(meanAmpSmall')'; % [ROIs x stimuli], amplitudes averaged across small pupil trials
meanAmpLarge = cellfun(@nanmean, ampLarge, repmat({2},size(amplitudes,1),1), ...
    'UniformOutput', false);
meanAmpLarge = cell2mat(meanAmpLarge')'; % [ROIs x stimuli], amplitudes averaged across small pupil trials
% mean amplitudes after shuffling stimulus labels
shuffleAmpSmall = cellfun(@nanmean, ampShuffleSmall, ...
    repmat({2},size(amplitudes,1),1), 'UniformOutput', false);
shuffleAmpSmall = cell2mat(shuffleAmpSmall); % [(stimuli*ROIs) x shuffles]
shuffleAmpSmall = reshape(shuffleAmpSmall, length(dirct), [], shuffles); % [stimuli x ROIs x shuffles] mean amplitudes after shuffling stimulus labels
shuffleAmpSmall = permute(shuffleAmpSmall, [2 1 3]); % [ROIs x stimuli x shuffles] 
shuffleAmpLarge = cellfun(@nanmean, ampShuffleLarge, ...
    repmat({2},size(amplitudes,1),1), 'UniformOutput', false);
shuffleAmpLarge = cell2mat(shuffleAmpLarge); % [(stimuli*ROIs) x shuffles]
shuffleAmpLarge = reshape(shuffleAmpLarge, length(dirct), [], shuffles); % [stimuli x ROIs x shuffles] mean amplitudes after shuffling stimulus labels
shuffleAmpLarge = permute(shuffleAmpLarge, [2 1 3]); % [ROIs x stimuli x shuffles] 
% inverte responses of suppressed ROIs
meanAmpSmall(isSuppr==1,:) = -meanAmpSmall(isSuppr==1,:);
meanAmpLarge(isSuppr==1,:) = -meanAmpLarge(isSuppr==1,:);
shuffleAmpSmall(isSuppr==1,:,:) = -shuffleAmpSmall(isSuppr==1,:,:);
shuffleAmpLarge(isSuppr==1,:,:) = -shuffleAmpLarge(isSuppr==1,:,:);
% set responses below baseline to zero
meanAmpSmall(meanAmpSmall<0) = 0;
meanAmpLarge(meanAmpLarge<0) = 0;
shuffleAmpSmall(shuffleAmpSmall<0) = 0;
shuffleAmpLarge(shuffleAmpSmall<0) = 0;
% standardize responses so they sum to one
meanAmpSmall = meanAmpSmall ./ nansum(meanAmpSmall,2);
meanAmpLarge = meanAmpLarge ./ nansum(meanAmpLarge,2);
shuffleAmpSmall = shuffleAmpSmall ./ nansum(shuffleAmpSmall,2);
shuffleAmpLarge = shuffleAmpLarge ./ nansum(shuffleAmpLarge,2);
% Determine DSIsSmall
vects = sum(dirVectors .* meanAmpSmall, 2);
shuffleVects = squeeze(sum(dirVectors .* shuffleAmpSmall, 2));
DSIsSmall = abs(vects);
nullDSIs = abs(shuffleVects);
p_DSISmall = sum(nullDSIs > DSIsSmall,2) ./ shuffles;
vects = sum(dirVectors .* meanAmpLarge, 2);
shuffleVects = squeeze(sum(dirVectors .* shuffleAmpLarge, 2));
DSIsLarge = abs(vects);
nullDSIs = abs(shuffleVects);
p_DSILarge = sum(nullDSIs > DSIsLarge,2) ./ shuffles;
% Determine OSIsSmall
vects = sum(oriVectors .* meanAmpSmall, 2);
shuffleVects = squeeze(sum(oriVectors .* shuffleAmpSmall, 2));
OSIsSmall = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSISmall = sum(nullOSIs > OSIsSmall,2) ./ shuffles;
vects = sum(oriVectors .* meanAmpLarge, 2);
shuffleVects = squeeze(sum(oriVectors .* shuffleAmpLarge, 2));
OSIsLarge = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSILarge = sum(nullOSIs > OSIsLarge,2) ./ shuffles;
% Determine pref directions
vects = squeeze(sum(dirVectors .* meanAmpSmall, 2)); % [neuron x small/large pupil]
prefDirsSmall = mod(angle(vects) ./ pi .* 180, 360);
vects = squeeze(sum(dirVectors .* meanAmpLarge, 2)); % [neuron x small/large pupil]
prefDirsLarge = mod(angle(vects) ./ pi .* 180, 360);
% Determine pref directions
vects = squeeze(sum(oriVectors .* meanAmpSmall, 2)); % [neuron x small/large pupil]
prefOrisSmall = mod(angle(vects) ./ pi .* 180, 360);
vects = squeeze(sum(oriVectors .* meanAmpLarge, 2)); % [neuron x small/large pupil]
prefOrisLarge = mod(angle(vects) ./ pi .* 180, 360);

%% Figure S2B (scatter: corr. with running vs corr. with pupil during gratings)
figure
plot(rhosRunGratings, rhosPupilGratings, 'k.', 'MarkerSize', 6)
ind = ~any(isnan([rhosRunGratings, rhosPupilGratings]),2);
[rho, pVal] = corr(rhosRunGratings(ind), rhosPupilGratings(ind));
xlabel('Correlation with running')
ylabel('Correlation with pupil')
title(sprintf('n = %d, rho = %.3f, p = %.2e', sum(ind), rho, pVal))
axis square equal
hold on
axis([-0.5 0.9 -0.5 0.9])
ax = gca;
ax.Box = 'off';

%% Figure S2C (proportion of large pupil trials for each grating direction)
lrgPpl = [];
drctn = [];
for k = 1:length(indDatasets)
    lp = largePupil{indDatasets(k)};
    lrgPpl = [lrgPpl; lp(:)];
    d = directions{indDatasets(k)};
    d(isnan(d)) = [];
    d = repmat(d, 1, size(lp,2));
    drctn = [drctn; d(:)];
end

p = anova1(lrgPpl, drctn, 'off');

bins = unique(drctn);
means = NaN(1, length(bins));
ses = NaN(1, length(bins));
for b = 1:length(bins)
    means(b) = nanmean(lrgPpl(drctn==bins(b)));
    ses(b) = nanstd(lrgPpl(drctn==bins(b))) ./ sqrt(sum(drctn==bins(b)));
end
means(end+1) = means(1);
ses(end+1) = ses(1);
bins(end+1) = 360;
figure
hold on
plot(bins, means, 'k', 'LineWidth', 2)
plot(bins, means+ses, 'k--')
plot(bins, means-ses, 'k--')
ylim([0.3 .5])
xlim([0 360])
xlabel('Direction of stimulus')
ylabel('Proportion of large pupil trials')
title(sprintf('ANOVA: p = %.3f', p))
set(gca, 'XTick', 0:90:360)

%% Figure S2D (scatter: OSIsSmall during small vs large pupil)
tbl = table(OSIsSmall, OSIsLarge, subjects, dataset, 'VariableNames', ...
    {'OSI_small', 'OSI_large', 'mouse', 'session'});
lme = fitlme(tbl, 'OSI_large ~ -1 + OSI_small + (-1 + OSI_small | session) + (-1 + OSI_small | mouse)');

edges = 0 : 1/100 : 1;
bins = edges(1:end-1)+diff(edges(1:2));
m = hot(200);
m = m(1:180,:);
N = histcounts2(OSIsSmall, OSIsLarge, edges, edges);
[xout, yout, zout] = prepareSurfaceData(bins, bins, N);
f = fit([xout, yout], zout, 'linearinterp');
densities = f([OSIsSmall, OSIsLarge]);
figure
hold on
scatter(OSIsSmall, OSIsLarge, [], densities, 'filled')
plot([0 1], [0 1], 'Color', [1 1 1].*0.8, 'LineWidth', 2)
plot([0 1], [0 1] .* fixedEffects(lme), 'r', 'LineWidth', 2)
colormap(m)
c = colorbar;
c.Label.String = 'density';
axis square
axis([0 1 0 1])
set(gca, 'XTick', [0 1], 'YTick', [0 1])
xlabel('OSI (small pupil)')
ylabel('OSI (large pupil)')
title(sprintf('n = %d, slope: %.3f, p = %.2e', ...
    sum(~any(isnan([OSIsSmall, OSIsLarge]),2)), fixedEffects(lme), coefTest(lme)))

%% Figure S2E (scatter: preferred directions during small vs large pupil)
ind = p_DSISmall < 0.05 & p_DSILarge < 0.05 & ~isnan(prefDirsSmall);
diffs = abs(prefDirsSmall(ind) - prefDirsLarge(ind));
diffs(diffs>180) = 360 - diffs(diffs>180);
figure
hold on
scatter(prefDirsSmall(ind), prefDirsLarge(ind), 'k', 'filled')
plot([0 360], [0 360], '--', 'Color', lines(1), 'LineWidth', 2)
axis square
axis([0 360 0 360])
set(gca, 'XTick', 0:90:360, 'YTick', 0:90:360)
xlabel('Preferred direction (small pupil)')
ylabel('Preferred direction (large pupil)')
title(sprintf('n = %d (%.1f%% with <20 degrees difference)', sum(ind), ...
    sum(diffs<20)/sum(ind)*100))

%% Figure S2F (scatter: preferred orientations during small vs large pupil)
ind = p_OSISmall < 0.05 & p_OSILarge < 0.05 & ~isnan(prefOrisSmall);
a = prefOrisSmall(ind)./180.*2.*pi - prefOrisLarge(ind)./180.*2.*pi;
j = a > pi;
a(j) = 2*pi - a(j);
j = a < - pi;
a(j) = 2*pi + a(j);
p = 1e-10;
while true
    h = circ_mtest(a, 0, p);
    if ~h
        break
    end
    p = p * 10;
end
diffs = abs(prefOrisSmall(ind) - prefOrisLarge(ind));
diffs(diffs>180) = 360 - diffs(diffs>180);
figure
hold on
scatter(prefOrisSmall(ind), prefOrisLarge(ind), 'k', 'filled')
plot([0 360], [0 360], '--', 'Color', lines(1), 'LineWidth', 2)
axis square
axis([0 360 0 360])
set(gca, 'XTick', 0:90:360, 'YTick', 0:90:360)
xlabel('Preferred orientation (small pupil)')
ylabel('Preferred orientation (large pupil)')
title(sprintf('n = %d (p < %.1e, circular paired t-test)', sum(ind), p))