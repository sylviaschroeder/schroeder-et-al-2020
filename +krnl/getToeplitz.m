function [A, numSamples, windowTimes]  = getToeplitz(time, eventTimes, ...
    evWindows, vectors, vecWindows, normalise)
% time          [t x 1]
% eventTimes    {1 x eventTypes}, each entry: [events x 1], times of each
%               event
% evWindows     {1 x eventTypes}, each entry: [start end] times of kernel
%               window relative to event onset
% vectors       {1 x continuousTypes}, each entry: [t x 1], values of
%               continuous measures
% vecWindows    {1 x continuousTypes}, each entry: [start end] times of
%               kernel window
% (normalise)   true or false (default); if true, STD of each predictor,
%               i.e. each column of toeplitz matrix is set to 1
%
% A             [t x predictors], toeplitz matrix
% numSamples    [1 x eventTypes+continuousTypes], length of each kernel
% windowTimes   {1 x eventTypes+continuousTypes}, kernel time relative to
%               event

if nargin < 6
    normalise = false;
end

Fs = 1/median(diff(time)); % frame rate

startOffset = zeros(1, length(evWindows));
numSamples = zeros(1, length(evWindows));
windowTimes = cell(1, length(evWindows));
for w = 1:length(evWindows)
    startOffset(w) = round(evWindows{w}(1)*Fs);
    numSamples(w) = round(diff(evWindows{w})*Fs);
    windowTimes{w} = (startOffset(w) + (0:numSamples(w)-1)) ./ Fs;
end

A = []; % the toeplitz matrix

for ev = 1:length(eventTimes)
    evSam = find(hist(eventTimes{ev}, time));
    eventFrames = evSam + startOffset(ev);
    firstColumn = zeros(length(time), 1);
    firstColumn(eventFrames(eventFrames>0)) = 1;
    firstRow = zeros(1, length(windowTimes{ev}));
    firstRow(1) = firstColumn(1);
    firstRow(1 - eventFrames(eventFrames<1 & ...
        -eventFrames<size(firstRow,2)-1)) = 1;
    A = [A, toeplitz(firstColumn, firstRow)];
end

if nargin > 3 && ~isempty(vectors)
    ns = zeros(1, length(vecWindows));
    wt = cell(1, length(vecWindows));
    for v = 1:length(vectors)
        win = round(vecWindows{v}(1)*Fs) : round(vecWindows{v}(2)*Fs);
        ns(v) = length(win);
        wt{v} = win ./ Fs;
        a = zeros(length(time), length(win));
        if ~isempty(vectors{v})
            for w = 1:length(win)
                a(max(1,1+win(w)) : min(length(time),length(time)+win(w)), end+1-w) = ...
                    vectors{v}(max(1,1-win(w)) : min(length(time),length(time)-win(w)));
            end
        end
        A = [A, a];
    end
    numSamples = [numSamples, ns];
    windowTimes = [windowTimes, wt];
end

if normalise
    s = std(A, 0, 1);
    s(s == 0) = 1;
    A = A ./ s;
end