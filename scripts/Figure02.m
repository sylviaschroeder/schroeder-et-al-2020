%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish';
folderTools = 'C:\STORAGE\workspaces';
folderThisRepo = 'C:\dev\workspace\schroeder-et-al-2020';

%% Parameters
% Preprocess neural data
smoothStd = 0.25; % in sec
sigma = 1; % in sec

%% Examples
examplesCorr = {'SS078', '2017-09-28', 1, 108; ...
            'SS078', '2017-09-28', 1, 215};

examplesTun = {'SS076', '2017-10-04', 1, 137; ...
            'SS077', '2017-10-05', 1, 24; ...
            'SS069', '2016-10-13', 3, 4; ...
            'SS077', '2017-10-03', 1, 56};

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Figure 2A,B (bouton traces during darkness and gratings)
planes = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pRois._ss_2pPlanes.npy'));
ids = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pRois.ids.npy'));
traces = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pCalcium.dff.npy'));
time = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pCalcium.timestamps.npy'));
delays = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_2pPlanes.delay.npy'));
running = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_running.speed.npy'));
runningTime = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_running.timestamps.npy'));
pupil = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\eye.diameter.npy'));
pupilTime = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\eye.timestamps.npy'));
gratingTimes = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_grating.intervals.npy'));
rhosDark = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_corrsRunning.rhosDark.npy'));
rhosGratings = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_corrsRunning.rhosGratings.npy'));
recTimeGratings = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_recordings.gratings_intervals.npy'));
recTimeDark = readNPY(fullfile(folderBase, 'boutons', examplesCorr{1,1}, ...
    examplesCorr{1,2}, '001\_ss_recordings.darkness_intervals.npy'));

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

r_p = mean([rhosDark rhosGratings],2);
[r_p_sorted,order] = sort(r_p,'descend');
ids_sorted = ids(order,:);
cols = 'rb';
stimuli = {'Gratings', 'Darkness'};
recTimes = [recTimeGratings; recTimeDark];
presTimes = [500 1200; 2330 3030];
for exp = 1:2
    startInd = find(time == recTimes(exp,1));
    endInd = find(time == recTimes(exp,2));
    tr = traces(startInd:endInd,:);
    filtered = filteredTraces(startInd:endInd,:);
    zscored = (filtered - nanmean(filtered,1)) ./ nanstd(filtered,0,1) ./ ...
        size(filtered,1).^0.5;
    zscored = zscored(:,order);
    
    r = runningNew(startInd:endInd);
    p = pupilNew(startInd:endInd);
    t = time(startInd:endInd);
    if ~all(isnan(p))
        ind = ~isnan(p);
        p = interp1(t(ind), p(ind), t, 'pchip');
    end
    
    cm = flip(gray);
    ax = [0 0 0];
    figure('Position',[15 45 1530 940])
    subplot(5,1,1)
    hold on
    if exp == 1 % gratings
        h = [0 0 0];
        h(1) = plot(t, p, 'Color', [0 0.7 0.5]);
        h(2) = plot(t, r./2-3, 'Color', [0 0 .7]);
        numStim = size(gratingTimes,1);
        stimT = [repmat(gratingTimes(:,1)',2,1);NaN(1,numStim)];
        s = [zeros(1,numStim); ones(1,numStim).*5; NaN(1,numStim)];
        h_ = plot(stimT,s-15,'k');
        h(3) = h_(1);
        leg = legend(h, 'Running','Pupil','Stimulus');
    else
        h = plot(t, r./2-3, 'Color', [0 0 .7]);
        leg = legend(h, 'Running');
    end
    title(sprintf('%s', stimuli{exp}))
    ylim([-20 40])
    set(gca,'box','off')
    ax(1) = gca;
    leg.Position = [0.92,0.83,0.05,0.09];
    subplot(5,1,2:4)
    imagesc(t([1 end]),[1 size(zscored,2)], zscored', ...
        prctile(zscored(:),[5 95]))
    colormap(cm)
    set(gca,'box','off')
    ax(2) = gca;
    subplot(5,1,5)
    hold on
    for ex = 1:size(examplesCorr,1)
        idx = planes == examplesCorr{ex,3} & ids == examplesCorr{ex,4};
        plot(t, smooth(traces(startInd:endInd,idx),5)-10*ex, cols(ex))
    end
    ax(3) = gca;
    linkaxes(ax,'x')
    xlim(t([find(t>presTimes(exp,1),1) find(t>presTimes(exp,2),1)]))
    xlabel('Time (s)')
    set(gca,'box','off')
end

%% Prepare for Figs. 2C-F
% Collect relevant variables from correlation data
subjects = {};
dates = {};
planes = [];
ids = [];
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
    end
end

pValsRunDark = sum(nullsRunDark < rhosRunDark, 2) ./ size(nullsRunDark, 2);
ind = pValsRunDark > 0.5;
pValsRunDark(ind) = 1 - pValsRunDark(ind);
pValsRunDark = 2 .* pValsRunDark; % two-sided test
pValsRunDark(pValsRunDark==0) = 1/size(nullsRunDark,2);

pValsRunGratings = sum(nullsRunGratings < rhosRunGratings, 2) ./ size(nullsRunGratings, 2);
ind = pValsRunGratings > 0.5;
pValsRunGratings(ind) = 1 - pValsRunGratings(ind);
pValsRunGratings = 2 .* pValsRunGratings; % two-sided test
pValsRunGratings(pValsRunGratings==0) = 1/size(nullsRunDark,2);

%% Figure 2C (correlations with running during darkness)
ind = isnan(rhosRunDark);
rhos = rhosRunDark;
rhos(ind) = [];
x = sort(rhos, 'ascend');
y = (1:length(x))' ./ length(x);
x = [-1; x; 1];
y = [0; y; 1];

nulls = nullsRunDark;
nulls(ind,:) = [];
xNull = sort(nulls, 1, 'ascend');
xNull = sort(xNull, 2, 'ascend');
limNull = prctile(xNull, [2.5 97.5], 2);
limNull = [[-1 -1]; limNull; [1 1]];
xNull = median(xNull, 2);
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [-1; xNull; 1];
yNull = [0; yNull; 1];

rhosEx = NaN(size(examplesCorr,1),1);
for ex = 1:size(examplesCorr,1)
    idx = strcmp(subjects, examplesCorr{ex,1}) & ...
        strcmp(dates, examplesCorr{ex,2}) & planes == examplesCorr{ex,3} & ...
        ids == examplesCorr{ex,4};
    rhosEx(ex) = rhosRunDark(idx);
end

cols = 'rb';
figure
hold on
h = [0 0];
h(2) = fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
    'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
h(1) = plot(x, y, 'k', 'LineWidth', 2);
heights = interp1(x, y, rhosEx);
for ex = 1:length(rhosEx)
    plot(rhosEx(ex), heights(ex), 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', cols(ex))
end
xlim([-0.52 0.52])
xlabel('Correlation with running')
ylabel('Proportion of boutons')
title(sprintf('Darkness (n = %d)', length(rhos)))
legend(h, {'boutons','shifted'}, 'Location', 'NorthWest')

%% Figure 2D (correlations with running during gratings)
ind = isnan(rhosRunGratings);
rhos = rhosRunGratings;
rhos(ind) = [];
x = sort(rhos, 'ascend');
y = (1:length(x))' ./ length(x);
x = [-1; x; 1];
y = [0; y; 1];

nulls = nullsRunGratings;
nulls(ind,:) = [];
xNull = sort(nulls, 1, 'ascend');
xNull = sort(xNull, 2, 'ascend');
limNull = prctile(xNull, [2.5 97.5], 2);
limNull = [[-1 -1]; limNull; [1 1]];
xNull = median(xNull, 2);
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [-1; xNull; 1];
yNull = [0; yNull; 1];

rhosEx = NaN(size(examplesCorr,1),1);
for ex = 1:size(examplesCorr,1)
    idx = strcmp(subjects, examplesCorr{ex,1}) & ...
        strcmp(dates, examplesCorr{ex,2}) & planes == examplesCorr{ex,3} & ...
        ids == examplesCorr{ex,4};
    rhosEx(ex) = rhosRunGratings(idx);
end

cols = 'rb';
figure
hold on
h = [0 0];
h(2) = fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
    'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
h(1) = plot(x, y, 'k', 'LineWidth', 2);
heights = interp1(x, y, rhosEx);
for ex = 1:length(rhosEx)
    plot(rhosEx(ex), heights(ex), 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', cols(ex))
end
xlim([-0.52 0.52])
xlabel('Correlation with running')
ylabel('Proportion of boutons')
title(sprintf('Gratings (n = %d)', length(rhos)))
legend(h, {'boutons','shifted'}, 'Location', 'NorthWest')

%% Figure 2E (scatterplot: correlations with running during darkness vs during gratings)
% Scatter plot
figure
hold on
plot(rhosRunDark, rhosRunGratings, 'k.', 'MarkerSize', 6)
cols = 'rb';
for ex = 1:size(examplesCorr,1)
    idx = strcmp(subjects, examplesCorr{ex,1}) & ...
        strcmp(dates, examplesCorr{ex,2}) & planes == examplesCorr{ex,3} & ...
        ids == examplesCorr{ex,4};
    plot(rhosRunDark(idx), rhosRunGratings(idx), [cols(ex) '.'], 'MarkerSize', 30);
end
ind = ~any(isnan([rhosRunDark, rhosRunGratings]), 2);
[rho,p] = corr(rhosRunDark(ind), rhosRunGratings(ind));
xlabel(sprintf('Correlation with running\nduring darkness'))
ylabel(sprintf('Correlation with running\nduring gratings'))
title(sprintf('n = %d, rho = %.3f, p = %.2e', sum(ind), rho, p))
axis([-0.65 0.85 -0.65 0.85])
ax = gca;
ax.Box = 'off';
axis square

edges = -0.65 : 0.1 : 0.85;
bins = edges(2:end)-0.05;
% Histogram for darkness
valid = ~isnan(rhosRunDark);
n1 = histcounts(rhosRunDark(pValsRunDark<0.05 & valid),edges)';
n2 = histcounts(rhosRunDark(pValsRunDark>=0.05 & valid),edges)';
figure;
b = bar(bins, [n1,n2], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
xlabel(sprintf('Correlation with running\nduring darkness'))
ylabel('#Boutons')
xlim(edges([1 end]))
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p \geq 0.05')
% Histogram for gratings
valid = ~isnan(rhosRunGratings);
n1 = histcounts(rhosRunGratings(pValsRunGratings<0.05 & valid),edges)';
n2 = histcounts(rhosRunGratings(pValsRunGratings>=0.05 & valid),edges)';
figure;
b = bar(bins, [n1,n2], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
xlabel(sprintf('Correlation with running\nduring gratings'))
ylabel('#Boutons')
xlim(edges([1 end]))
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p \geq 0.05')

%% Figure 2F (correlations with pupil during gratings)
ind = isnan(rhosPupilGratings);
rhos = rhosPupilGratings;
rhos(ind) = [];
x = sort(rhos, 'ascend');
y = (1:length(x))' ./ length(x);
x = [-1; x; 1];
y = [0; y; 1];

nulls = nullsPupilGratings;
nulls(ind,:) = [];
xNull = sort(nulls, 1, 'ascend');
xNull = sort(xNull, 2, 'ascend');
limNull = prctile(xNull, [2.5 97.5], 2);
limNull = [[-1 -1]; limNull; [1 1]];
xNull = median(xNull, 2);
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [-1; xNull; 1];
yNull = [0; yNull; 1];

rhosEx = NaN(size(examplesCorr,1),1);
for ex = 1:size(examplesCorr,1)
    idx = strcmp(subjects, examplesCorr{ex,1}) & ...
        strcmp(dates, examplesCorr{ex,2}) & planes == examplesCorr{ex,3} & ...
        ids == examplesCorr{ex,4};
    rhosEx(ex) = rhosPupilGratings(idx);
end

cols = 'rb';
figure
hold on
h = [0 0];
h(2) = fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
    'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
h(1) = plot(x, y, 'k', 'LineWidth', 2);
heights = interp1(x, y, rhosEx);
for ex = 1:length(rhosEx)
    plot(rhosEx(ex), heights(ex), 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', cols(ex))
end
xlim([-0.52 0.52])
xlabel('Correlation with pupil')
ylabel('Proportion of boutons')
title(sprintf('Gratings (n = %d)', length(rhos)))
legend(h, {'boutons','shifted'}, 'Location', 'NorthWest')

%% Figure 2G (example tuning curves, small and large pupil)
cols = 'kr';
for ex = 1:size(examplesTun,1)
    planes = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pRois._ss_2pPlanes.npy'));
    ids = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_2pRois.ids.npy'));
    directions = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingID.directions.npy'));
    largePupil = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_gratingTrials.largePupil.npy'));
    curvesSmall = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_tuning.curvesSmall.npy'));
    curvesLarge = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
        examplesTun{ex,2}, '001\_ss_tuning.curvesLarge.npy'));
    curves = cat(3, curvesSmall, curvesLarge);
    amplitudes = readNPY(fullfile(folderBase, 'boutons', examplesTun{ex,1}, ...
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

%% Prepare for Figs. 2H,I
% Collect relevant variables from tuning data
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
    end
end

% invert sign of responses of suppressed cells
minima(isSuppr==1,:) = -minima(isSuppr==1,:);
maxima(isSuppr==1,:) = -maxima(isSuppr==1,:);
means(isSuppr==1,:) = -means(isSuppr==1,:);
nullMinima(isSuppr==1,:,:) = -nullMinima(isSuppr==1,:,:);
nullMaxima(isSuppr==1,:,:) = -nullMaxima(isSuppr==1,:,:);
nullMeans(isSuppr==1,:,:) = -nullMeans(isSuppr==1,:,:);

shuffles = 1000;
directions = 0:30:330;
dirVectors = exp(directions./180.*pi .* 1i);
oriVectors = exp(directions./180.*2.*pi .* 1i);

ampSmall = amplitudes; % {nROIs x 1}, each entry: [nStim x nReps]
ampLarge = amplitudes;
for n = 1:size(amplitudes,1)
    ampSmall{n}(largePupil{n}) = NaN;
    ampLarge{n}(~largePupil{n}) = NaN;
end
meanAmpSmall = cellfun(@nanmean, ampSmall, repmat({2},size(amplitudes,1),1), ...
    'UniformOutput', false);
meanAmpSmall = cell2mat(meanAmpSmall')'; % [ROIs x stimuli], amplitudes averaged across small pupil trials
meanAmpLarge = cellfun(@nanmean, ampLarge, repmat({2},size(amplitudes,1),1), ...
    'UniformOutput', false);
meanAmpLarge = cell2mat(meanAmpLarge')'; % [ROIs x stimuli], amplitudes averaged across small pupil trials

% inverte responses of suppressed ROIs
meanAmpSmall(isSuppr==1,:) = -meanAmpSmall(isSuppr==1,:);
meanAmpLarge(isSuppr==1,:) = -meanAmpLarge(isSuppr==1,:);
% set responses below baseline to zero
meanAmpSmall(meanAmpSmall<0) = 0;
meanAmpLarge(meanAmpLarge<0) = 0;
% standardize responses so they sum to one
meanAmpSmall = meanAmpSmall ./ nansum(meanAmpSmall,2);
meanAmpLarge = meanAmpLarge ./ nansum(meanAmpLarge,2);

% Determine DSIs
vects = sum(dirVectors .* meanAmpSmall, 2);
DSIsSmall = abs(vects);

vects = sum(dirVectors .* meanAmpLarge, 2);
DSIsLarge = abs(vects);

%% Figure 2H (response modulations)
cols = lines(4);
binSize = 20;
mini = -80;
maxi = 80;
modFun = @(a,b) (b-a)./((abs(a)+abs(b)) ./ 2) .* 100;

mx = maxima;
ind = all(isnan(maxima),2);
mx(ind,:) = means(ind,:);
nmx = nullMaxima;
nmx(ind,:,:) = nullMeans(ind,:,:);
respMod = modFun(mx(:,1), mx(:,2));
nullMod = modFun(squeeze(nmx(:,1,:)), squeeze(nmx(:,2,:)));
confInt = prctile(nullMod, [2.5 97.5], 2);
sgnfcnt = respMod < confInt(:,1) | respMod > confInt(:,2);

% Histogram
figure
bins = mini:binSize:maxi;
edges = [bins-binSize/2, maxi+binSize/2];
n1 = histcounts(respMod(sgnfcnt), edges);
n2 = histcounts(respMod(~sgnfcnt), edges);
b = bar(bins, [n1',n2'], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
hold on
plot(nanmean(respMod(sgnfcnt & respMod<0)), 1100, 'vk', 'MarkerFaceColor', 'k')
plot(nanmean(respMod(sgnfcnt & respMod>0)), 1100, 'vk', 'MarkerFaceColor', 'k')
xlim(edges([1 end]))
title(sprintf('n = %d', sum(~isnan(respMod))))
xlabel('Response modulation (%)')
ylabel('#Boutons')
legend(b, 'p < 0.05', 'p \geq 0.05')
ax = gca;
ax.Box = 'off';
ax.XTick = [mini 0 maxi];

% cumulative distribution
figure
hold on
h = [0 0];
plot([0 0], [0 1], 'k')
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
ylabel('Proportion of boutons')
title(sprintf('n = %d', sum(~isnan(respMod))))
legend(h, {'boutons','shifted'}, 'Location', 'NorthWest')
legend('boxoff')
set(gca, 'XTick', [mini 0 maxi]);

%% Figure 2I (DSIs small vs large pupil)
[~,~,sbj] = unique(subjects);
[~,~,dt] = unique(dates);
[~,~,dataset] = unique([sbj, dt], 'rows');
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