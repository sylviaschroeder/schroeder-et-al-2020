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
examplesRF = {'SS096', '2018-03-08', 70; ...
            'SS098', '2018-03-16', 65};

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Figure 3A (receptive fields)
titles = {'ON field','OFF field'};
for ex = 1:size(examplesRF,1)
    ids = readmatrix(fullfile(folderBase, 'opticTract', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\clusters.uuids.csv'));
    rfs = readNPY(fullfile(folderBase, 'opticTract', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_rf.maps.npy'));
    pos = readNPY(fullfile(folderBase, 'opticTract', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\_ss_rfDescr.edges.npy'));
    
    unit = ids == examplesRF{ex,3};
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
for ex = 1:size(examplesRF,1)
    ids = readmatrix(fullfile(folderBase, 'opticTract', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\clusters.uuids.csv'));
    % NEED TO GET ALL WAVEFORMS! not just mean
    wfs = readNPY(fullfile(folderBase, 'opticTract', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\clusters.waveforms.npy'));
    coord = readNPY(fullfile(folderBase, 'opticTract', examplesRF{ex,1}, ...
        examplesRF{ex,2}, '001\channels.localCoordinates.npy'));
    
    unit = ids == examplesRF{ex,3};
    % find channel with largest spike
    [~, peakChan] = max(max(abs(mean(squeeze(wfs(unit,:,:)),1)),[],3),[],2);
    plotChans = max(1,peakChan-numChans) : min(length(sp(probe).xcoords),peakChan+numChans);
    ycoords = unique(sp(probe).ycoords(plotChans));
    mini = min(ycoords);
    maxi = max(ycoords);
    chanDist = median(diff(ycoords));
    
    figure('Position', [1145 42 754 1074]);
    hold on
    h = [0 0];
    for c = 1:2
        H = arrayfun(@(x)plot(sp(probe).xcoords(x)+0.3*(1:size(wf.waveFormsMean,3))', ...
            sp(probe).ycoords(x)+.5.*squeeze(wf.waveFormsMean(c,x,:)), cols(c)), plotChans);
        h(c) = H(1);
    end
    ylim([mini - 0.5*chanDist, maxi + 1.5 * chanDist])
    ylabel('Depth (um)')
    set(gca, 'XTick', [])
    legend(h, {'stationary','running'}, 'Location', 'NorthEast')
    title(sprintf('%s %s %s, unit %d', db(dataset).subject, ...
        db(dataset).date, sp(probe).name, unit))
end