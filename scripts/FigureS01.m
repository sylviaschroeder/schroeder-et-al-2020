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

%% Prepare for Figs. S1C-F
% Collect relevant variables from tuning data
subjects = {};
dates = {};
planes = [];
ids = [];
prefDirs = [];
isSuppr = [];
amplitudes = {};
largePupil = {};
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
        prefDirs = [prefDirs; parsS(:,1)];
        isSuppr = [isSuppr; readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_tuning.isSuppressed.npy'))];
        amp = readNPY(fullfile(folderBase, 'boutons', name, ...
            date, '001\_ss_gratingTrials.amplitudes.npy'));
        amplitudes = [amplitudes; permute(mat2cell(amp, ...
            size(amp,1), size(amp,2), ones(1,n)), [3 1 2])];
        largePupil = [largePupil; repmat({readNPY(fullfile(folderBase, ...
            'boutons', name, date, ...
            '001\_ss_gratingTrials.largePupil.npy'))}, n, 1)];
        
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
[~,~,sbj] = unique(subjects);
[~,~,dt] = unique(dates);
[~,~,dataset] = unique([sbj, dt], 'rows');

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

% Determine DSIs and OSIs
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

%% Figure S1C (Preferred direction for each dataset)
binSize = 22.5;
minNum = 10;
edges = 0:binSize:360;
bins = edges(1:end-1) + binSize/2;
dSets = unique(dataset);
hists = NaN(max(dSets), length(bins));
for d = 1:length(dSets)
    ind = dataset==dSets(d) & (p_DSI<0.05 | p_OSI<0.05);
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
ylabel('Proportion of boutons')

%% Figure S1D (ON/OFF index for boutons driven and suppressed by gratings)
% Note: due to bug in code used for paper, there are small discrapencies
% between paper figure and figure here (results don't change qualitatively)
groups = validRF & [isSuppr == -1, isSuppr == 1];
groupNames = {sprintf('driven (%d)', sum(groups(:,1) & ~isnan(OnOffRatios))), ...
    sprintf('suppressed (%d)', sum(groups(:,2) & ~isnan(OnOffRatios)))};

tbl = table(OnOffRatios(validRF), nominal(isSuppr(validRF)), ...
    dataset(validRF), subjects(validRF), ...
    'VariableNames', {'OnOff','isSuppr','session','mouse'});
lme = fitlme(tbl, 'OnOff ~ isSuppr + (1|session) + (isSuppr-1|session)', ...
    'DummyVarCoding','effects');

figure
hold
plot([0 0], [0 1], 'k')
h = zeros(1, size(groups,2));
colors = [0 0 0; 0.5 0.5 0.5];
b = fixedEffects(lme);
means = [sum(b) b(1)-b(2)];
for g = 1:size(groups,2)
    x = sort(OnOffRatios(groups(:,g) & ~isnan(OnOffRatios)), 'ascend');
    y = (1:sum(groups(:,g) & ~isnan(OnOffRatios))) ./ ...
        sum(groups(:,g) & ~isnan(OnOffRatios));
    x = [-1; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(means(g), 1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
end
ylim([0 1.01])
legend(h, groupNames, 'Location', 'SouthEast')
xlabel('ON/OFF index')
ylabel('Proportion of boutons')
title(sprintf('p = %.4f', coefTest(lme)))

%% Figure S1E (DSI for ON, OFF and ON+OFF boutons)
groups = [validRF & OnOffRatios>onThr, ...
    validRF & OnOffRatios<offThr, ...
    validRF & OnOffRatios<=onThr & OnOffRatios>=offThr];
groupNames = {'ON', 'OFF', 'ON+OFF'};

RFtypes = cell(length(validRF),1);
for g = 1:3
    [RFtypes{groups(:,g)}] = deal(groupNames{g});
end
% check whether ON, OFF, ON+OFF units have different mean DSI
tbl = table(DSIs(validRF), RFtypes(validRF), dataset(validRF), ...
    subjects(validRF), 'VariableNames', {'DSI','RF','session','mouse'});
lme = fitlme(tbl, 'DSI ~ RF + (RF|session)', 'DummyVarCoding', 'effects');

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
    ind = groups(:,g) & ~isnan(DSIs);
    x = sort(DSIs(ind), 'ascend');
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
xlim([0 0.8])
ylim([0 1.01])
legend(h, lbls, 'Location', 'SouthEast')
xlabel('DSI')
ylabel('Proportion of boutons')
title(sprintf('ANOVA: p = %.2e', anova(lme).pValue(2)))
annotation('textbox', [0.5 0.3 0.4 0.2], 'String', text, ...
    'LineStyle', 'none')

%% Figure S1F (OSI for ON, OFF and ON+OFF boutons)
groups = [validRF & OnOffRatios>onThr, ...
    validRF & OnOffRatios<offThr, ...
    validRF & OnOffRatios<=onThr & OnOffRatios>=offThr];
groupNames = {'ON', 'OFF', 'ON+OFF'};

RFtypes = cell(length(validRF),1);
for g = 1:3
    [RFtypes{groups(:,g)}] = deal(groupNames{g});
end
% check whether ON, OFF, ON+OFF units have different mean OSI
tbl = table(OSIs(validRF), RFtypes(validRF), dataset(validRF), ...
    subjects(validRF), 'VariableNames', {'OSI','RF','session','mouse'});
lme = fitlme(tbl, 'OSI ~ RF + (RF|session)', 'DummyVarCoding', 'effects');

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
    ind = groups(:,g) & ~isnan(OSIs);
    x = sort(OSIs(ind), 'ascend');
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
xlim([0 0.6])
ylim([0 1.01])
legend(h, lbls, 'Location', 'SouthEast')
xlabel('OSI')
ylabel('Proportion of boutons')
title(sprintf('ANOVA: p = %.2e', anova(lme).pValue(2)))
annotation('textbox', [0.5 0.3 0.4 0.2], 'String', text, ...
    'LineStyle', 'none')