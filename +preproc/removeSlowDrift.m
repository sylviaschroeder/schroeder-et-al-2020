function [filteredTraces, smoothed] = removeSlowDrift(traces, ...
    samplingRate, window, percentile)

if nargin < 3 || isempty(window)
    window = 100; % in s
end
if nargin < 4 || isempty(percentile)
    percentile = 5;
end

smoothed = zeros(size(traces));
n = round(window * samplingRate);
if mod(n,2) == 0
    n = n+1;
end
nBefore = floor((n-1)/2);
nAfter = n - nBefore - 1;

a = gcp('nocreate');
if isempty(a)
    for k = 1:size(traces,1)
        tmpTraces = traces(max(1,k-nBefore) : min(size(traces,1),k+nAfter),:);
        smoothed(k,:) = prctile(tmpTraces, percentile);
    end
else
    parfor k = 1:size(traces,1)
        tmpTraces = traces(max(1,k-nBefore) : min(size(traces,1),k+nAfter),:);
        smoothed(k,:) = prctile(tmpTraces, percentile);
    end
end
% smoothed = sgolayfilt(smoothed, 3, n);

filteredTraces = traces - smoothed;