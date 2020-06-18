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

% running speed and firing rates
binSizeRun = 1/7.5;
sigma = 1;
sig = round(sigma / binSizeRun);
win = normpdf(-5*sig : 5*sig, 0, sig);
highPassWindow = 180; % in sec; to high-pass filter firing rates
prctileFilter = 8;
runThreshold = 1;
minRunTime = 3;

%% Examples
examples = {'SS096', '2018-03-08', 70; ...
            'SS098', '2018-03-16', 65};

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Figure 3A (receptive fields)
titles = {'ON field','OFF field'};
for ex = 1:size(examples,1)
    ids = readmatrix(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\clusters.uuids.csv'));
    rfs = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_rf.maps.npy'));
    pos = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_rfDescr.edges.npy'));
    
    unit = ids == examples{ex,3};
    rf = squeeze(rfs(unit,:,:,:,:));
    rf(:,:,:,2) = -rf(:,:,:,2);
    squW = diff(pos(1:2)) / size(rf,2);
    squH = diff(pos(3:4)) / size(rf,1);
    [mx,mxTime] = max(max(reshape(permute(abs(rf),[1 2 4 3]),[],size(rf,3)),[],1));
    figure
    for f = 1:2
        subplot(2,1,f)
        imagesc([pos(1)+squW/2 pos(2)-squW/2], [pos(3)+squH/2 pos(4)-squH/2], ...
            rf(:,:,mxTime,f),[-mx mx])
        axis image
        set(gca, 'box', 'off', 'XTick', pos(1:2), 'YTick', [pos(3) 0 pos(4)], ...
            'YTickLabel', [-pos(3) 0 -pos(4)])
        colormap(gca, colormaps{f})
        ylabel(titles{f})
        colorbar
        if f == 1
            title(sprintf('RGC axon %d', ex))
        end
    end
end

%% Figure 3B (waveforms)
numChans = 10;
cols = 'kr';
for ex = 1:size(examples,1)
    ids = readmatrix(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\clusters.uuids.csv'));
    wfs = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\clusters.waveforms.npy'));
    coord = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\channels.localCoordinates.npy'));
    
    unit = ids == examples{ex,3};
    % find channel with largest spike
    [~, peakChan] = max(max(abs(mean(squeeze(wfs(unit,:,:,:)),3)),[],2),[],1);
    plotChans = max(1,peakChan-numChans) : min(size(coord,1),peakChan+numChans);
    ycoords = unique(coord(plotChans,2));
    mini = min(ycoords);
    maxi = max(ycoords);
    chanDist = median(diff(ycoords));
    
    figure('Position', [1145 42 754 1074]);
    hold on
    h = [0 0];
    for c = 1:2
        H = arrayfun(@(x)plot(coord(x,1)+0.3*(1:size(wfs,3))', ...
            coord(x,2)+.5.*squeeze(wfs(unit,x,:,c)), cols(c)), plotChans);
        h(c) = H(1);
    end
    ylim([mini - 0.5*chanDist, maxi + 1.5 * chanDist])
    ylabel('Depth (um)')
    set(gca, 'XTick', [])
    legend(h, {'stationary','running'}, 'Location', 'NorthEast')
    title(sprintf('RGC axon %d', ex))
end

%% Figure 3C (running speed and firing rate during darkness)
for ex = 1:size(examples,1)
    ids = readmatrix(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\clusters.uuids.csv'));
    spikeTimes = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\spike.times.npy'));
    clusters = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\spike.clusters.npy'));
    validIntervals = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_validTimes.intervals.npy'));
    validClusters = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_validTimes.clusters.npy'));
    darkTime = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_darkness.intervals.npy'));
    runSpeed = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_running.speed.npy'));
    runTime = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_running.timestamps.npy'));
    
    binEdges = darkTime(1) : binSizeRun : darkTime(2);
    timeBins = binEdges(1:end-1) + binSizeRun/2;
    vInt = validIntervals(validClusters == unit,:);
    
    unit = examples{ex,3};
    st = spikeTimes(clusters == unit);
    st(st<darkTime(1) | st>darkTime(2)) = [];
    [sr,~,b] = histcounts(st, binEdges);
    sr = [ones(1,length(win)) .* mean(sr(1:min(length(sr),length(win)))), ...
        sr, ones(1,length(win)) .* ...
        mean(sr(end-min(length(sr),length(win))+1:end))] ./ binSizeRun;
    sr = conv(sr, win, 'same');
    sr = sr(length(win)+1 : end-length(win));
    [sr, smoothed] = preproc.removeSlowDrift(sr', 1/binSizeRun, ...
        highPassWindow, prctileFilter);
    sr = sr' + mean(smoothed);
    
    run = interp1(runTime, runSpeed, timeBins, 'pchip'); % run speed during darkness only
    run = [ones(1,length(win)) .* mean(run(1:min(length(run),length(win)))), ...
        run, ones(1,length(win)) .* mean(run(end-min(length(run),length(win))+1:end))];
    run = conv(run, win, 'same');
    run = run(length(win)+1 : end-length(win));
    isRunning = run > runThreshold;
    
    sta = find(diff(isRunning)==1);
    sto = find(diff(isRunning)==-1);
    if sta(1)>sto(1)
        sta = [1, sta];
    end
    if sta(end)>sto(end)
        sto(end+1) = length(beh);
    end
    ind = find(sta(2:end) - sto(1:end-1) < minRunTime/binSizeRun);
    sta(ind+1) = [];
    sto(ind) = [];
    ind = (sto - sta) >= minRunTime/binSizeRun;
    sta = sta(ind);
    sto = sto(ind);

    if ~isempty(vInt)
        sr2 = NaN(size(sr));
        run2 = NaN(size(run));
        for k = 1:size(vInt,1)
            ind = timeBins>=vInt(k,1) & timeBins<=vInt(k,2);
            sr2(ind) = sr(ind);
            run2(ind) = run(ind);
        end
        sr = sr2;
        run = run2;
    end
    
    figure('Position', [4 678 1914 420]);
    subplot(2,1,1)
    hold on
    subplot(2,1,2)
    hold on
    for k = 1:length(sta)
        subplot(2,1,1)
        fill(timeBins([sta(k) sto(k) sto(k) sta(k)]), [[1 1].*min(run), [1 1].*max(run)], 'k', ...
            'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.1)
        subplot(2,1,2)
        fill(timeBins([sta(k) sto(k) sto(k) sta(k)]), [[1 1].*min(sr), [1 1].*max(sr)], 'k', ...
            'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.1)
    end
    ax = [0 0];
    subplot(2,1,1)
    plot(timeBins, run, 'k')
    ylabel('Running speed (cm/s)')
    title(sprintf('RGC axon %d', ex))
    axis tight
    ax(1) = gca;
    subplot(2,1,2)
    plot(timeBins, sr, 'k')
    xlabel('Time (s)')
    ylabel('Firing rate (sp/s)')
    axis tight
    ax(2) = gca;
    linkaxes(ax, 'x')
    xlim(timeBins([1 end]))
end

%% Figure 3D (cross-correlograms)
maxLag = 30;
for ex = 1:size(examples,1)
    ids = readmatrix(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\clusters.uuids.csv'));
    crossCorrs = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_crossCorrs.values.npy'));
    lags = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_crossCorrs.timestamps.npy'));
    nullCorrs = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_crossCorrs.nullValues.npy'));
    
    unit = ids == examples{ex,3};
    null = prctile(nullCorrs(:,:,unit), [2.5 97.5], 2);
    real = crossCorrs(:,unit);
    mini = floor(min([null(:);real(:)]) * 10) / 10;
    maxi = ceil(max([null(:);real(:)]) * 10) / 10;
    figure
    hold on
    fill([lags; flip(lags)], [null(:,1); flip(null(:,2))], 'k', ...
        'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
    plot(lags, real, 'k', 'LineWidth', 1)
    plot(lags([1 end]), [0 0], 'k:')
    plot([0 0], [mini maxi], 'k:')
    xlim([-maxLag maxLag])
    ylim([mini maxi])
    set(gca, 'YTick', mini:0.1:maxi, 'XTick', [-maxLag 0 maxLag])
    xlabel('Time lag (s)')
    ylabel('Cross-correlation')
    title(sprintf('RGC axon %d', ex))
end

%% Figure 3E (scatter: correlation vs p-value)
pVals = [];
rhos =[];

subjDirs = dir(fullfile(folderBase, 'opticTract', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'opticTract', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        crossCorrs = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_crossCorrs.values.npy'));
        lags = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_crossCorrs.timestamps.npy'));
        nullCorrs = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_crossCorrs.nullValues.npy'));
    
        [~,zeroInd] = min(abs(lags));
        for n = 1:size(crossCorrs,2)
            rhos(end+1,1) = crossCorrs(zeroInd,n);
            null = nullCorrs(zeroInd,:,n);
            p = sum(null > rhos(end)) / length(null);
            if p > 0.5
                p = 1-p;
            end
            pVals(end+1,1) = 2 * p;
        end
    end
end
pVals(pVals==0) = 1e-4;

figure
scatter(rhos, pVals, 'k', 'filled')
hold on
plot([-.4 0.4],[1 1].*0.05, 'k')
set(gca, 'YScale', 'log', 'YDir', 'reverse', 'YMinorTick', false, 'XTick', [-.4 0 .4])
xlabel('Correlation with running')
ylabel('p-value')
title(['n = ' num2str(length(rhos))])