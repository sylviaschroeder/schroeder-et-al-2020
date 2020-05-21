%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish';
folderTools = 'C:\STORAGE\workspaces';
folderThisRepo = 'C:\dev\workspace\schroeder-et-al-2020';

%% Parameters
RFtypes = {'ON','OFF','ON+OFF'};

% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

% to group into ON, OFF and ON+OFF
onThr = 0.5; % ON field is at least three times stronger than OFF field
offThr = -0.5; % OFF field is at least three times stronger than ON field

% colormaps
red = [1 0 .5];
% green = [0 1 .5];
blue = [0 .5 1];
% orange = [1 .5 0];
black = [1 1 1].*0.5;
grad = linspace(0,1,100)';
reds = red.*flip(grad) + [1 1 1].*grad;
% greens = green.*flip(grad) + [1 1 1].*grad;
blacks = black.*flip(grad) + [1 1 1].*grad;
% cm_ON = [greens; flip(reds(1:end-1,:),1)];
cm_ON = [blacks; flip(reds(1:end-1,:),1)];
blues = blue.*flip(grad) + [1 1 1].*grad;
% oranges = orange.*flip(grad) + [1 1 1].*grad;
% cm_OFF = [oranges; flip(blues(1:end-1,:),1)];
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

%% Figure 1F
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

%% Figure 1G
% Collect relevant variables from visual noise data
subjects = {};
dates = {};
planes = [];
ids = [];
rfs = [];
evTotal = [];
evStim = [];
evRun = [];
lambdasStim = [];
pValues = [];
OnOffValues = [];

subjDirs = dir(fullfile(folderBase, 'boutons', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'boutons', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        folder = fullfile(folderBase, 'boutons', name, date, '001');
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
        rfs = [rfs; rfs_exp];
        evTotal = [evTotal; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_rf.explVars.npy'))];
        evStim = [evStim; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_rf.explVarsStim.npy'))];
        evRun = [evRun; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_rf.explVarsRunning.npy'))];
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

%% Figure 1H & I
buffer = 2; % in sec (before and after stim period)
cols = [0 0 0; 1 0 0];
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


%% Figure 1J