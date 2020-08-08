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
examplesRF = {'SS041', '2015-04-23', 4, 58, -1; ...
            'SS048', '2015-12-02', 3, 62, 1; ...
            'SS044', '2015-04-28', 3, 188, 1; ...
            'SS041', '2015-04-23', 4, 43, -1; ...
            'SS044', '2015-05-29', 2, 72, -1; ...
            'SS044', '2015-04-28', 2, 204, 1};

examplesTun = {'SS038', '2015-02-17', 1, 127; ...
            'SS041', '2015-04-23', 2, 38; ...
            'SS041', '2015-04-23', 2, 151; ...
            'SS038', '2015-02-17', 1, 203};

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(genpath(fullfile(folderTools, 'CircStat2012a')))
addpath(fullfile(folderThisRepo))

%% Figure S4A (receptive fields of example neurons)
titles = {'ON field','OFF field'};
for ex = 1:size(examplesRF,1)
    roiPlanes = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_2pRois._ss_2pPlanes.npy'));
    ids = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_2pRois.ids.npy'));
    rfs = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_rf.maps.npy'));
    pos = readNPY(fullfile(folderBase, 'sc neurons 2p', examplesRF{ex,1}, ...
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

%% Prepare for Figures S4B-V
subjects = {};
dates = {};
planes = [];
ids = [];

evStim = [];
lambdasStim = [];
pValues = [];
OnOffValues = [];

numSh = 200;
prefDirs = [];
minima = [];
maxima = [];
means = [];
nullMinima = [];
nullMaxima = [];
nullMeans = [];
isGad = [];
isSuppr = [];
amplitudes = {};
largePupil = {};
directions = {};

rhosRunGray = [];
nullsRunGray = [];
rhosPupilGray = [];
nullsPupilGray = [];
rhosPupilGratings = [];
nullsPupilGratings = [];

subjDirs = dir(fullfile(folderBase, 'sc neurons 2p', 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    fprintf('Subject %s (%d of %d)\n', name, subj, length(subjDirs))
    dateDirs = dir(fullfile(folderBase, 'sc neurons 2p', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        
        p = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_2pRois._ss_2pPlanes.npy'));
        n = length(p);
        planes = [planes; p];
        ids = [ids; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
            date, '001\_ss_2pRois.ids.npy'))];
        subjects = [subjects; repmat({name}, n, 1)];
        dates = [dates; repmat({date}, n, 1)];
        
        % correlations
        if isfile(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsRunning.rhosGrayScreen.npy'))
            rho = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsRunning.rhosGrayScreen.npy'));
            null = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsRunning.nullRhosGrayScreen.npy'));
        else
            rho = NaN(n,1);
            null = NaN(n,500);
        end
        rhosRunGray = [rhosRunGray; rho];
        nullsRunGray = [nullsRunGray; null];
        if isfile(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsPupil.rhosGrayScreen.npy'))
            rho = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsPupil.rhosGrayScreen.npy'));
            null = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsPupil.nullRhosGrayScreen.npy'));
        else
            rho = NaN(n,1);
            null = NaN(n,500);
        end
        rhosPupilGray = [rhosPupilGray; rho];
        nullsPupilGray = [nullsPupilGray; null];
        if isfile(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsPupil.rhosGratings.npy'))
            rho = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsPupil.rhosGratings.npy'));
            null = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_corrsPupil.nullRhosGratings.npy'));
        else
            rho = NaN(n,1);
            null = NaN(n,500);
        end
        rhosPupilGratings = [rhosPupilGratings; rho];
        nullsPupilGratings = [nullsPupilGratings; null];
        
        % tuning data
        mi = NaN(n,2);
        ma = NaN(n,2);
        mn = NaN(n,2);
        nmi = NaN(n,2,numSh);
        nma = NaN(n,2,numSh);
        nmn = NaN(n,2,numSh);
        if isfile(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_tuning.parametersSmall.npy'))
            parsS = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_tuning.parametersSmall.npy'));
            parsL = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_tuning.parametersLarge.npy'));
            pars = cat(3, parsS, parsL);
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
            pd = parsS(:,1);
            pd(isnan(parsS(:,2))) = NaN;
            prefDirs = [prefDirs; pd];
            isGad = [isGad; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_2pRois.isGad.npy'))];
            isSuppr = [isSuppr; readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_tuning.isSuppressed.npy'))];
            amp = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_gratingTrials.amplitudes.npy'));
            amplitudes = [amplitudes; permute(mat2cell(amp, ...
                size(amp,1), size(amp,2), ones(1,n)), [3 1 2])];
            largePupil = [largePupil; repmat({readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_gratingTrials.largePupil.npy'))}, n, 1)];
            directions = [directions; repmat({readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_gratingID.directions.npy'))}, n, 1)];
            
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
        else
            prefDirs = [prefDirs; NaN(n,1)];
            isGad = [isGad; NaN(n,1)];
            isSuppr = [isSuppr; NaN(n,1)];
            amplitudes = [amplitudes; cell(n,1)];
            largePupil = [largePupil; cell(n,1)];
            directions = [directions; cell(n,1)];
        end
        minima = [minima; mi];
        maxima = [maxima; ma];
        means = [means; mn];
        nullMinima = [nullMinima; nmi];
        nullMaxima = [nullMaxima; nma];
        nullMeans = [nullMeans; nmn];
        
        % receptive fields
        oov = NaN(n,2);
        if ~isfile(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_rf.maps.npy'))
            ev = NaN(n,1);
            lam = NaN(n,1);
            pv = NaN(n,1);
        else
            rfs = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_rf.maps.npy'));
            ev = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_rf.explVarsStim.npy'));
            lam = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_rf.lambdasStim.npy'));
            pv = readNPY(fullfile(folderBase, 'sc neurons 2p', name, ...
                date, '001\_ss_rf.pValues.npy'));
            
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

% determine response modulation and tuning depth modulation
modFun = @(a,b) (b-a)./((abs(a)+abs(b))./2).*100;
mx = maxima;
nmx = nullMaxima;
ind = all(isnan(mx),2);
mx(ind,:) = means(ind,:);
nmx(ind,:,:) = nullMeans(ind,:,:);
respMod = modFun(mx(:,1), mx(:,2));
nullRespMod = modFun(squeeze(nmx(:,1,:)), squeeze(nmx(:,2,:)));
depth = abs(maxima - minima);
nullDepth = abs(nullMaxima - nullMinima);
depthMod = modFun(depth(:,1), depth(:,2));
nullDepthMod = modFun(squeeze(nullDepth(:,1,:)), squeeze(nullDepth(:,2,:)));

% determine whether neurons are tuned
isTuned = ~any(isnan(maxima),2);

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
p_DSISmall(isnan(DSIsSmall)) = NaN;
vects = sum(dirVectors .* meanAmpLarge, 2);
shuffleVects = squeeze(sum(dirVectors .* shuffleAmpLarge, 2));
DSIsLarge = abs(vects);
nullDSIs = abs(shuffleVects);
p_DSILarge = sum(nullDSIs > DSIsLarge,2) ./ shuffles;
p_DSILarge(isnan(DSIsLarge)) = NaN;
% Determine OSIsSmall
vects = sum(oriVectors .* meanAmpSmall, 2);
shuffleVects = squeeze(sum(oriVectors .* shuffleAmpSmall, 2));
OSIsSmall = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSISmall = sum(nullOSIs > OSIsSmall,2) ./ shuffles;
p_OSISmall(isnan(OSIsSmall)) = NaN;
vects = sum(oriVectors .* meanAmpLarge, 2);
shuffleVects = squeeze(sum(oriVectors .* shuffleAmpLarge, 2));
OSIsLarge = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSILarge = sum(nullOSIs > OSIsLarge,2) ./ shuffles;
p_OSILarge(isnan(OSIsLarge)) = NaN;
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

%% Figure S4B (histogram of ON/OFF indices)
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
ylabel('#Neurons')
title(sprintf('n = %d', sum(validRF)))

%% Figure S4C (scatter of DSIs vs OSIs)
% Test significance for separation between high OSIs and high DSIs
goodDSIs = DSIsSmall(~isnan(DSIsSmall) & isTuned & (p_DSISmall<.05|p_OSISmall<.05));
goodOSIs = OSIsSmall(~isnan(OSIsSmall) & isTuned & (p_DSISmall<.05|p_OSISmall<.05));
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
mx = ceil(max([DSIsSmall; OSIsSmall]) / binSize) * binSize;
colors = [1 0 0; 0 0 1; .5 0 .5];
inds = [p_DSISmall<.05&p_OSISmall>=.05, p_DSISmall>=.05&p_OSISmall<.05, ...
    p_DSISmall<.05&p_OSISmall<.05];
h = [0 0 0];
figure
hold on
plot([0 mx], [0 mx], 'k')
for k = 3:-1:1
    h(k) = scatter(DSIsSmall(isTuned & inds(:,k)), ...
        OSIsSmall(isTuned & inds(:,k)), [], colors(k,:), 'filled');
end
for ex = 1:size(examplesTun,1)
    idx = strcmp(subjects, examplesTun{ex,1}) & ...
        strcmp(dates, examplesTun{ex,2}) & planes==examplesTun{ex,3} & ...
        ids==examplesTun{ex,4};
    if p_DSISmall(idx)>=0.05 && p_OSISmall(idx)>=0.05
        continue
    end
    scatter(DSIsSmall(idx), OSIsSmall(idx), [], cols(ex,:), 'LineWidth', 2)
end
legend(h,'dir-sel','ori-sel','dir&ori')
axis square
axis([0 mx 0 mx])
set(gca, 'box','off','XTick',0:.1:mx,'YTick',0:.1:mx)
xlabel('DSI')
ylabel('OSI')
title(sprintf('n = %d (DSI and OSI are different, p = %.2e)', ...
    sum(~isnan(DSIsSmall) & isTuned & (p_DSISmall<0.05|p_OSISmall<0.05)), pDiff))

%% Figure S4D (ON/OFF index for excitatory and inhibitory neurons)
% Note: figure in paper was (unnecessarily) restricted to neurons that were
% responsive to gratings; here all neurons are included (results don't
% change qualitatively)
groups = validRF & [isGad == -1, isGad == 1];
groupNames = {sprintf('excitatory (%d)', sum(groups(:,1) & ~isnan(OnOffRatios))), ...
    sprintf('inhibitory (%d)', sum(groups(:,2) & ~isnan(OnOffRatios)))};

tbl = table(OnOffRatios(any(groups,2)), nominal(isGad(any(groups,2))), ...
    dataset(any(groups,2)), subjects(any(groups,2)), ...
    'VariableNames', {'OnOff','inhibitory','session','mouse'});
lme = fitlme(tbl, 'OnOff ~ inhibitory + (1|session) + (1|mouse)', ...
    'DummyVarCoding','effects');

figure
hold
plot([0 0], [0 1], 'k')
h = zeros(1, size(groups,2));
colors = [0 0 0; 0.5 0.5 0.5];
b = fixedEffects(lme);
m = [sum(b) b(1)-b(2)];
for g = 1:size(groups,2)
    x = sort(OnOffRatios(groups(:,g) & ~isnan(OnOffRatios)), 'ascend');
    y = (1:sum(groups(:,g) & ~isnan(OnOffRatios))) ./ ...
        sum(groups(:,g) & ~isnan(OnOffRatios));
    x = [-1; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(m(g), 1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
end
ylim([0 1.01])
legend(h, groupNames, 'Location', 'SouthEast')
xlabel('ON/OFF index')
ylabel('Proportion of neurons')
title(sprintf('p = %.4f', coefTest(lme)))

%% Figure S4E (DSI for excitatory and inhibitory neurons)
groups = [isGad == -1, isGad == 1];
groupNames = {sprintf('excitatory (%d)', sum(groups(:,1) & ~isnan(DSIsSmall))), ...
    sprintf('inhibitory (%d)', sum(groups(:,2) & ~isnan(DSIsSmall)))};

tbl = table(DSIsSmall(any(groups,2)), nominal(isGad(any(groups,2))), ...
    dataset(any(groups,2)), subjects(any(groups,2)), ...
    'VariableNames', {'DSI','inhibitory','session','mouse'});
lme = fitlme(tbl, 'DSI ~ inhibitory + (1|session) + (inhibitory-1|session)', ...
    'DummyVarCoding','effects');

figure
hold
plot([0 0], [0 1], 'k')
h = zeros(1, size(groups,2));
colors = [0 0 0; 0.5 0.5 0.5];
b = fixedEffects(lme);
m = [sum(b) b(1)-b(2)];
for g = 1:size(groups,2)
    x = sort(DSIsSmall(groups(:,g) & ~isnan(DSIsSmall)), 'ascend');
    y = (1:sum(groups(:,g) & ~isnan(DSIsSmall))) ./ ...
        sum(groups(:,g) & ~isnan(DSIsSmall));
    x = [-1; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(m(g), 1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
end
xlim([0 max(DSIsSmall)])
ylim([0 1.01])
legend(h, groupNames, 'Location', 'SouthEast')
xlabel('DSI')
ylabel('Proportion of neurons')
title(sprintf('p = %.4f', coefTest(lme)))

%% Figure S4F (OSI for excitatory and inhibitory neurons)
groups = [isGad == -1, isGad == 1];
groupNames = {sprintf('excitatory (%d)', sum(groups(:,1) & ~isnan(OSIsSmall))), ...
    sprintf('inhibitory (%d)', sum(groups(:,2) & ~isnan(OSIsSmall)))};

tbl = table(OSIsSmall(any(groups,2)), nominal(isGad(any(groups,2))), ...
    dataset(any(groups,2)), subjects(any(groups,2)), ...
    'VariableNames', {'DSI','inhibitory','session','mouse'});
lme = fitlme(tbl, 'DSI ~ inhibitory + (1|session) + (inhibitory-1|session) + (1|mouse)', ...
    'DummyVarCoding','effects');

figure
hold
plot([0 0], [0 1], 'k')
h = zeros(1, size(groups,2));
colors = [0 0 0; 0.5 0.5 0.5];
b = fixedEffects(lme);
m = [sum(b) b(1)-b(2)];
for g = 1:size(groups,2)
    x = sort(OSIsSmall(groups(:,g) & ~isnan(OSIsSmall)), 'ascend');
    y = (1:sum(groups(:,g) & ~isnan(OSIsSmall))) ./ ...
        sum(groups(:,g) & ~isnan(OSIsSmall));
    x = [-1; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(m(g), 1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
end
xlim([0 max(OSIsSmall)])
ylim([0 1.01])
legend(h, groupNames, 'Location', 'SouthEast')
xlabel('OSI')
ylabel('Proportion of neurons')
title(sprintf('p = %.4f', coefTest(lme)))

%% Figure S4G (ON/OFF index for neurons driven- and suppressed-by-gratings)
groups = validRF & [isSuppr == -1, isSuppr == 1];
groupNames = {sprintf('driven (%d)', sum(groups(:,1) & ~isnan(OnOffRatios))), ...
    sprintf('suppressed (%d)', sum(groups(:,2) & ~isnan(OnOffRatios)))};

tbl = table(OnOffRatios(any(groups,2)), nominal(isSuppr(any(groups,2))), ...
    dataset(any(groups,2)), subjects(any(groups,2)), ...
    'VariableNames', {'OnOff','isSuppr','session','mouse'});
lme = fitlme(tbl, 'OnOff ~ isSuppr + (1|session) + (isSuppr-1|session)', ...
    'DummyVarCoding','effects');

figure
hold
plot([0 0], [0 1], 'k')
h = zeros(1, size(groups,2));
colors = [0 0 0; 0.5 0.5 0.5];
b = fixedEffects(lme);
m = [sum(b) b(1)-b(2)];
for g = 1:size(groups,2)
    x = sort(OnOffRatios(groups(:,g) & ~isnan(OnOffRatios)), 'ascend');
    y = (1:sum(groups(:,g) & ~isnan(OnOffRatios))) ./ ...
        sum(groups(:,g) & ~isnan(OnOffRatios));
    x = [-1; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(m(g), 1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
end
ylim([0 1.01])
legend(h, groupNames, 'Location', 'SouthEast')
xlabel('ON/OFF index')
ylabel('Proportion of neurons')
title(sprintf('p = %.4f', coefTest(lme)))

%% Figure S4H (DSI for ON, OFF and ON+OFF neurons)
groups = [validRF & OnOffRatios>onThr, ...
    validRF & OnOffRatios<offThr, ...
    validRF & OnOffRatios<=onThr & OnOffRatios>=offThr];
groupNames = {'ON', 'OFF', 'ON+OFF'};

RFtypes = cell(length(validRF),1);
for g = 1:3
    [RFtypes{groups(:,g)}] = deal(groupNames{g});
end
% check whether ON, OFF, ON+OFF units have different mean DSI
tbl = table(DSIsSmall(validRF), RFtypes(validRF), dataset(validRF), ...
    subjects(validRF), 'VariableNames', {'DSI','RF','session','mouse'});
lme = fitlme(tbl, 'DSI ~ RF + (RF-1|mouse)', 'DummyVarCoding', 'effects');

[b, bnames] = fixedEffects(lme);
bnames = bnames.Variables;
means = NaN(1, size(groups,2));
figure;
hold
h = zeros(1, size(groups,2));
lbls = cell(1, size(groups,2));
colors = [1 0 1; 0 1 1; 0.5 0.5 1];
text = '';
for g = 1:size(groups,2)
    ind = groups(:,g) & ~isnan(DSIsSmall);
    x = sort(DSIsSmall(ind), 'ascend');
    y = (1:sum(ind)) ./ sum(ind);
    x = [0; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    v = find(strcmp(bnames, ['RF_' groupNames{g}]));
    if ~isempty(v)
        means(g) = b(1) + b(v);
        if v == 2
            text = [text, sprintf('%s vs %s: p = %.2e\n', groupNames{g}, ...
                bnames{3}(4:end), coefTest(lme, [0 1 -1]))];
        end
    else
        means(g) = b(1) - sum(b(2:end));
        for g2 = g+1:size(groups,2)
            v2 = find(strcmp(bnames, ['RF_' groupNames{g2}]));
            inds = zeros(1,3);
            inds(v2) = 1;
            text = [text, sprintf('%s vs %s: p = %.2e\n', groupNames{g}, ...
                groupNames{g2}, coefTest(lme, inds))];
        end
    end
    plot(means(g), 1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
    lbls{g} = sprintf('%s (n=%d)', groupNames{g}, sum(ind));
end
xlim([0 max(DSIsSmall)])
ylim([0 1.01])
legend(h, lbls, 'Location', 'SouthEast')
xlabel('DSI')
ylabel('Proportion of neurons')
title(sprintf('ANOVA: p = %.2e', anova(lme).pValue(2)))
annotation('textbox', [0.5 0.3 0.4 0.2], 'String', text, ...
    'LineStyle', 'none')

%% Figure S4I (OSI for ON, OFF and ON+OFF neurons)
groups = [validRF & OnOffRatios>onThr, ...
    validRF & OnOffRatios<offThr, ...
    validRF & OnOffRatios<=onThr & OnOffRatios>=offThr];
groupNames = {'ON', 'OFF', 'ON+OFF'};

RFtypes = cell(length(validRF),1);
for g = 1:3
    [RFtypes{groups(:,g)}] = deal(groupNames{g});
end
% check whether ON, OFF, ON+OFF units have different mean OSI
tbl = table(OSIsSmall(validRF), RFtypes(validRF), dataset(validRF), ...
    subjects(validRF), 'VariableNames', {'OSI','RF','session','mouse'});
lme = fitlme(tbl, 'OSI ~ RF', 'DummyVarCoding', 'effects');

[b, bnames] = fixedEffects(lme);
bnames = bnames.Variables;
means = NaN(1, size(groups,2));
figure;
hold
h = zeros(1, size(groups,2));
lbls = cell(1, size(groups,2));
colors = [1 0 1; 0 1 1; 0.5 0.5 1];
text = '';
for g = 1:size(groups,2)
    ind = groups(:,g) & ~isnan(OSIsSmall);
    x = sort(OSIsSmall(ind), 'ascend');
    y = (1:sum(ind)) ./ sum(ind);
    x = [0; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    v = find(strcmp(bnames, ['RF_' groupNames{g}]));
    if ~isempty(v)
        means(g) = b(1) + b(v);
        if v == 2
            text = [text, sprintf('%s vs %s: p = %.2e\n', groupNames{g}, ...
                bnames{3}(4:end), coefTest(lme, [0 1 -1]))];
        end
    else
        means(g) = b(1) - sum(b(2:end));
        for g2 = g+1:size(groups,2)
            v2 = find(strcmp(bnames, ['RF_' groupNames{g2}]));
            inds = zeros(1,3);
            inds(v2) = 1;
            text = [text, sprintf('%s vs %s: p = %.2e\n', groupNames{g}, ...
                groupNames{g2}, coefTest(lme, inds))];
        end
    end
    plot(means(g), 1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
    lbls{g} = sprintf('%s (n=%d)', groupNames{g}, sum(ind));
end
xlim([0 max(OSIsSmall)])
ylim([0 1.01])
legend(h, lbls, 'Location', 'SouthEast')
xlabel('OSI')
ylabel('Proportion of neurons')
title(sprintf('ANOVA: p = %.2e', anova(lme).pValue(2)))
annotation('textbox', [0.5 0.3 0.4 0.2], 'String', text, ...
    'LineStyle', 'none')

%% Figure S4J (Preferred direction for each dataset)
binSize = 22.5;
minNum = 10;
edges = 0:binSize:360;
bins = edges(1:end-1) + binSize/2;
dSets = unique(dataset);
hists = NaN(max(dSets), length(bins));
for d = 1:length(dSets)
    ind = dataset==dSets(d) & (p_DSISmall<0.05 | p_OSISmall<0.05);
    if sum(ind) < minNum
        continue
    end
    hists(dSets(d),:) = histcounts(prefDirs(ind), edges);
end
hists = hists ./ sum(hists, 2);
m = nanmean(hists,1);
se = nanstd(hists,0,1)./sqrt(sum(~isnan(hists),1));
figure
hold on
set(gca, 'ColorOrder', jet(max(dSets)))
plot(bins, hists, 'LineWidth', 2)
fill([bins flip(bins)], [m+se, flip(m-se)], 'k', 'EdgeColor', 'none', ...
    'FaceColor', 'k', 'FaceAlpha', 0.2)
plot(bins, m, 'k', 'LineWidth', 2)
xlim([0 360])
set(gca, 'XTick', 0:90:360, 'box', 'off')
xlabel('Preferred direction')
ylabel('Proportion of neurons')

%% Figure S4K (proportion of large pupil trials for each grating direction)
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
m = NaN(1, length(bins));
ses = NaN(1, length(bins));
for b = 1:length(bins)
    m(b) = nanmean(lrgPpl(drctn==bins(b)));
    ses(b) = nanstd(lrgPpl(drctn==bins(b))) ./ sqrt(sum(drctn==bins(b)));
end
m(end+1) = m(1);
ses(end+1) = ses(1);
bins(end+1) = 360;
figure
hold on
plot(bins, m, 'k', 'LineWidth', 2)
plot(bins, m+ses, 'k--')
plot(bins, m-ses, 'k--')
ylim([0.3 .507])
xlim([0 360])
xlabel('Direction of stimulus')
ylabel('Proportion of large pupil trials')
title(sprintf('ANOVA: p = %.3f', p))
set(gca, 'XTick', 0:90:360)

%% Figure S4L (scatterplot: correlations with pupil during gray screen vs during gratings)
figure
hold on
plot(rhosPupilGray, rhosPupilGratings, 'k.', 'MarkerSize', 6)
ind = ~any(isnan([rhosPupilGray, rhosPupilGratings]), 2);
[rho,p] = corr(rhosPupilGray(ind), rhosPupilGratings(ind));
xlabel(sprintf('Correlation with pupil\n(grey screen)'))
ylabel(sprintf('Correlation with pupil\n(visually driven)'))
title(sprintf('n = %d, rho = %.3f, p = %.2e', sum(ind), rho, p))
axis([-0.65 0.85 -0.65 0.85])
ax = gca;
ax.Box = 'off';
axis square

%% Figure S4M (scatterplot: correlations with pupil vs running (during gray screen)
figure
hold on
plot(rhosRunGray, rhosPupilGray, 'k.', 'MarkerSize', 6)
ind = ~any(isnan([rhosRunGray, rhosPupilGray]), 2);
[rho,p] = corr(rhosRunGray(ind), rhosPupilGray(ind));
xlabel(sprintf('Correlation with running\n(grey screen)'))
ylabel(sprintf('Correlation with pupil\n(grey screen)'))
title(sprintf('n = %d, rho = %.3f, p = %.2e', sum(ind), rho, p))
axis([-0.65 0.85 -0.65 0.85])
ax = gca;
ax.Box = 'off';
axis square

%% Figure S4N (DSIs small vs large pupil)
tbl = table(DSIsSmall, DSIsLarge, subjects, dataset, 'VariableNames', ...
    {'DSI_small', 'DSI_large', 'mouse', 'session'});
lme = fitlme(tbl, 'DSI_large ~ -1 + DSI_small + (-1 + DSI_small | session) + (-1 + DSI_small | mouse)');

edges = 0 : 1/100 : 1;
bins = edges(1:end-1)+diff(edges(1:2));
m = hot(200);
m = m(1:180,:);
N = histcounts2(DSIsSmall, DSIsLarge, edges, edges);
[xout, yout, zout] = prepareSurfaceData(bins, bins, N);
f = fit([xout, yout], zout, 'linearinterp');
densities = f([DSIsSmall, DSIsLarge]);
figure
hold on
scatter(DSIsSmall, DSIsLarge, [], densities, 'filled')
plot([0 1], [0 1], 'Color', [1 1 1].*0.8, 'LineWidth', 2)
plot([0 1], [0 1] .* fixedEffects(lme), 'r', 'LineWidth', 2)
colormap(m)
c = colorbar;
c.Label.String = 'density';
axis square
axis([0 1 0 1])
set(gca, 'XTick', [0 1], 'YTick', [0 1])
xlabel('DSI (small pupil)')
ylabel('DSI (large pupil)')
title(sprintf('n = %d, slope: %.3f, p = %.2e', ...
    sum(~any(isnan([DSIsSmall, DSIsLarge]),2)), fixedEffects(lme), coefTest(lme)))

%% Figure S4O (OSIs small vs large pupil)
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

%% Figure S4P (scatter: preferred directions during small vs large pupil)
ind = p_DSISmall < 0.05 & p_DSILarge < 0.05 & ~isnan(prefDirsSmall);
a = prefDirsSmall(ind)./180.*2.*pi - prefDirsLarge(ind)./180.*2.*pi;
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
title(sprintf('n = %d (mean diff.s = 0, p < %.1e, circular paired t-test)', sum(ind), p))

%% Figure S4Q (scatter: preferred orientations during small vs large pupil)
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
title(sprintf('n = %d (mean diff.s = 0, p < %.1e, circular paired t-test)', sum(ind), p))

%% Figure S4R (histogram and cum. distr. of tuning depth modulations)
cols = lines(4);
binSize = 20;
mini = -140;
maxi = 140;

confInt = prctile(nullDepthMod, [2.5 97.5], 2);
sgnfcnt = depthMod < confInt(:,1) | depthMod > confInt(:,2);

% Histogram
figure
bins = mini:binSize:maxi;
edges = [bins-binSize/2, maxi+binSize/2];
n1 = histcounts(depthMod(sgnfcnt), edges);
n2 = histcounts(depthMod(~sgnfcnt), edges);
b = bar(bins, [n1',n2'], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
xlim(edges([1 end]))
title(sprintf('n = %d', sum(~isnan(depthMod))))
xlabel('Tuning depth modulation (%)')
ylabel('#Neurons')
legend(b, 'p < 0.05', 'p \geq 0.05')
ax = gca;
ax.Box = 'off';
ax.XTick = [mini 0 maxi];

% cumulative distribution
figure
hold on
h = [0 0];
plot([0 0], [0 1], 'k')
ind = ~isnan(depthMod);
x = sort(depthMod(ind), 'ascend');
y = (1:sum(ind)) ./ sum(ind);
x = [-200; x; 200];
y = [0 y 1];
pseudo = nullDepthMod;
pseudo(~ind,:) = [];
xNull = sort(pseudo, 1, 'ascend');
xNull = sort(xNull, 2, 'ascend');
limNull = prctile(xNull, [2.5 97.5], 2);
limNull = [[-1 -1]; limNull; [1 1]];
yNull = (1:length(x)) ./ length(x);

h(2) = fill([limNull(:,1);flip(limNull(:,2))], [yNull, flip(yNull)], ...
    'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
h(1) = plot(x, y, 'k', 'LineWidth', 2);
xlim(edges([1 end]))
xlabel('Tuning depth modulation (%)')
ylabel('Proportion of neurons')
title(sprintf('n = %d', sum(~isnan(depthMod))))
legend(h, {'neurons','permuted'}, 'Location', 'NorthWest')
legend('boxoff')
set(gca, 'XTick', [mini 0 maxi]);

%% Figure S4S-V (distributions of correlations and tuning modulations for different cell types)
measures = [rhosPupilGray, rhosPupilGratings, respMod, depthMod];
nullMeasures = {nullsPupilGray, nullsPupilGratings, nullRespMod, nullDepthMod};
measureNames = {'Correlation with pupil during gray screens', ...
    'Correlation with pupil during gratings', ...
    'Response modulation (%)','Tuning depth modulation (%)'};
cellTypes = {{OnOffRatios>onThr & validRF; ...
    OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
    OnOffRatios<offThr & validRF}, ... % ON, ON+OFF, OFF
    {isGad==-1; isGad==1&isSuppr==-1; isGad==1&isSuppr==1}, ... %exc., inh.&driven, inh.&suppressed
    {any([p_OSISmall, p_OSILarge]<.05,2) & all([p_DSISmall, p_DSILarge]>=0.05,2); ...
    any([p_DSISmall, p_DSILarge]<.05,2) & all([p_OSISmall, p_OSILarge]>=0.05,2); ...
    any([p_DSISmall, p_DSILarge]<.05,2) & any([p_OSISmall, p_OSILarge]<.05,2)}}; % DS, OS, DS & OS
typeNames = {{'ON','ON+OFF','OFF'},{'exc','inh & driv','inh & sup'}, ...
    {'OS only','DS only','OS & DS'}};
colors = {[1 0 1; .5 .5 1; 0 1 1], [0.5 0 0; 1 .5 0; .75 1 0], [.5 1 0; 0 .5 1; .25 .75 .5]};
extremes = [-1 1; -1 1; -200 200; -200 200];
xLimits = [-0.5 0.5; -0.5 0.5; -60 60; -80 80];
valid = [repmat(~any(isnan([rhosPupilGray, rhosPupilGratings]), 2), 1, 2), ...
    ~isnan(respMod), ~isnan(depthMod)];
tests = {[2,6,3], [1,4,2], [7,5,8], [7,9,2]};

for m = 1:4
    for t = 1:3
        figure
        hold on
        plot([0 0], [0 1], 'k')
        h = zeros(1, length(cellTypes{t}));
        lbls = cell(1, length(cellTypes{t}));
        types = NaN(size(valid,1),1);
        for type = 1:length(cellTypes{t})
            ind = cellTypes{t}{type} & valid(:,m);
            hs = plots.plotCumHist(gca, measures(ind,m), ...
                nullMeasures{m}(ind,:), extremes(m,:), colors{t}(type,:));
            plot(mean(measures(ind,m)), 1.02, 'v', 'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', colors{t}(type,:))
            h(type) = hs(1);
            lbls{type} = sprintf('%s (n=%d)', typeNames{t}{type}, sum(ind));
            types(cellTypes{t}{type}) = type;
        end
        xlim(xLimits(m,:))
        ylim([0 1.05])
        legend(h, lbls, 'Location', 'SouthEast')
        xlabel(measureNames{m})
        ylabel('Proportion of neurons')
        
        tbl = table(measures(valid(:,m),m), categorical(types(valid(:,m))), ...
            subjects(valid(:,m)), categorical(dataset(valid(:,m))), ...
            'VariableNames', {'measure', 'type', 'mouse', 'session'});
        switch tests{m}(t)
            case 1
                lme = fitlme(tbl, 'measure ~ type + (type|session)', 'DummyVarCoding', 'effects');
            case 2
                lme = fitlme(tbl, 'measure ~ type + (1|session) + (type-1|session)', 'DummyVarCoding', 'effects');
            case 3
                lme = fitlme(tbl, 'measure ~ type + (1|session) + (1|mouse)', 'DummyVarCoding', 'effects');
            case 4
                lme = fitlme(tbl, 'measure ~ type + (1|session) + (type-1|mouse)', 'DummyVarCoding', 'effects');
            case 5
                lme = fitlme(tbl, 'measure ~ type + (1|session) + (1|mouse) + (type-1|mouse)', 'DummyVarCoding', 'effects');
            case 6
                lme = fitlme(tbl, 'measure ~ type + (1|session) + (type-1|session) + (type|mouse)', 'DummyVarCoding', 'effects');
            case 7
                lme = fitlme(tbl, 'measure ~ type + (1|session) + (type-1|session) + (1|mouse)', 'DummyVarCoding', 'effects');
            case 8
                lme = fitlme(tbl, 'measure ~ type + (type-1|session) + (1|mouse)', 'DummyVarCoding', 'effects');
            case 9
                lme = fitlme(tbl, 'measure ~ type + (type-1|session) + (1|mouse) + (type-1|mouse)', 'DummyVarCoding', 'effects');
        end
        stats = anova(lme);
        title(sprintf('p = %.4f (ANOVA)', stats.pValue(2)))
        if stats.pValue(2) < 0.05
            if length(cellTypes{t}) < 3
                fprintf('Pairwise p-vlaues for %s between %s, %s:\n', ...
                    measureNames{m}, typeNames{t}{:})
                fprintf('%.2e\n', coefTest(lme,[0 1]))
            else
                fprintf('Pairwise p-vlaues for %s between %s, %s, %s:\n', ...
                    measureNames{m}, typeNames{t}{:})
                fprintf('  %.2e\n', coefTest(lme,[0 1 0]), ...
                    coefTest(lme,[0 0 1]), coefTest(lme,[0 1 -1]))
            end
        end
    end
end