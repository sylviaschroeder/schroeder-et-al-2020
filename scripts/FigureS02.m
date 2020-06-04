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
addpath(fullfile(folderThisRepo))

%% Figure S2A
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

subjDirs = dir(fullfile(folderBase, 'boutons', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
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
            directions = [directions; repmat({readNPY(fullfile(folder, ...
                '_ss_gratingID.directions.npy'))}, n, 1)];
            
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
    end
end