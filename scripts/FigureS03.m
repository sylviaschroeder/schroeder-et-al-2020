%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish';
folderTools = 'C:\STORAGE\workspaces';
folderThisRepo = 'C:\dev\workspace\schroeder-et-al-2020';

%% Parameters
% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

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
examples = {'SS096', '2018-03-08', 70; ...
            'SS098', '2018-03-16', 65};

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(genpath(fullfile(folderTools, 'spikes')))
addpath(fullfile(folderThisRepo))

%% Figure S3B+D (responses to flickering monitor)
binSize = 0.001;
limits = [-.1 .4];
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
    flColor = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_flicker.color.npy'));
    flFreq = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_flicker.frequencies.npy'));
    flRep = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_flicker.repetition.npy'));
    flTimes = readNPY(fullfile(folderBase, 'opticTract', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_flicker.times.npy'));
    
    freqs = unique(flFreq);
    
    unit = examples{ex,3};
    reps = max(flRep);
    bins = limits(1):binSize:limits(2);
    bins = bins(1:end-1);
    
    flPSTH = NaN(reps,length(bins));
    for r = 1:reps
        ind = find(flFreq == freqs(end) & flRep==r); % flickers of fastest freq and repetiion r
        ind = ind(1:floor(length(ind)/2)*2); % only consider white/black flicker pairs (not single whites at end)
        ind = ind(1:2:end); % only whites
        flPSTH(r,:) = psthAndBA(spikeTimes(clusters==unit), ...
            flTimes(ind), limits, binSize);
    end
    
    onToOn = median(diff(flTimes(flFreq == freqs(end) & flColor == 1)));
    flicks = (-1:2)'.*onToOn;
    flicks = [flicks, flicks + onToOn/2];
    stim = [ones(size(flicks,1),1), zeros(size(flicks,1),1)];
    ax = [0 0];
    figure
    subplot(6,1,1)
    stairs(reshape(flicks',[],1), reshape(stim',[],1), 'k')
    ylim([-.1 1.1])
    ax(1) = gca;
    ylabel('Luminance')
    title(sprintf('Flickering monitor (RGC axon %d)', ex))
    
    subplot(6,1,2:6)
    stairs(bins, mean(flPSTH), 'k')
    ax(2) = gca;
    linkaxes(ax, 'x')
    xlim([-1/3 2] .* onToOn)
    set(ax(1), 'box', 'off')
    set(ax(2), 'box', 'off')
    xlabel('Time (s)')
    ylabel('Firing rate (sp/s)')
end

%% Figure S3E
subjDirs = dir(fullfile(folderBase, 'opticTract', 'SS*'));
rfNames = {'ON','OFF'};
maxLag = 30;
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, 'opticTract', name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        rfs = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_rf.maps.npy'));
        pos = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_rfDescr.edges.npy'));
        time = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_rfDescr.timestamps.npy'));
        crossCorrs = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_crossCorrs.values.npy'));
        lags = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_crossCorrs.timestamps.npy'));
        nullCorrs = readNPY(fullfile(folderBase, 'opticTract', name, ...
            date, '001\_ss_crossCorrs.nullValues.npy'));
        
        for unit = 1:size(rfs,1)
            rf = squeeze(rfs(unit,:,:,:,:));
            rf(:,:,:,2) = -rf(:,:,:,2);
            mx = max(abs(rf(:)));
            squW = diff(pos(1:2)) / size(rf,2);
            squH = diff(pos(3:4)) / size(rf,1);
            figure('Position', [10 560 1900 420])
            for f = 1:2
                for t = 1:2
                    subplot(1,5,(f-1)*2+3-t)
                    imagesc([pos(1)+squW/2 pos(2)-squW/2], [pos(3)+squH/2 pos(4)-squH/2], ...
                        rf(:,:,t,f),[-mx mx])
                    axis image
                    if f == 1 && t == 2
                        set(gca, 'box', 'off', 'XTick', [ceil(pos(1)) floor(pos(2))], ...
                            'YTick', [pos(3) 0 pos(4)], ...
                            'YTickLabel', [-pos(3) 0 -pos(4)])
                    else
                        set(gca, 'box', 'off', 'XTick', [], 'YTick', [])
                    end
                    colormap(gca, colormaps{f})
                    title(sprintf('%s (%.2fs)', rfNames{f}, -time(t)))
                    if t == 2
%                         if f == 1
%                             p = [0.259 0.34 0.008 0.348];
%                         else
%                             p = [0.584 0.34 0.008 0.348];
%                         end
%                         colorbar('Position', p)
                        if f == 1
                            p = [0.13 0.257 0.1237 0.018];
                        else
                            p = [0.4558 0.257 0.1237 0.018];
                        end
                        colorbar('Location', 'Southoutside', 'Position', p)
                    end
                end
            end
            
            null = prctile(nullCorrs(:,:,unit), [2.5 97.5], 2);
            real = crossCorrs(:,unit);
            mini = floor(min([null(:);real(:)]) * 10) / 10;
            maxi = ceil(max([null(:);real(:)]) * 10) / 10;
            subplot(1,5,5)
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
            set(gca, 'Position', [.79 0.25 0.124 0.5])
        end
    end
end