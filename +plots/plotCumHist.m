function h = plotCumHist(ax, data, null, xLimits, colors, alpha)

if nargin < 4
    xLimits = [-1 1];
end
if nargin < 5
    colors = [0 0 0; 0 0 0];
end
if size(colors,1) == 1
    colors = repmat(colors,2,1);
end
if nargin < 6
    alpha = 0.2;
end

axes(ax)
hold on
h = 0;
leg = {'data'};

if ~isempty(null)
    h = [0 0];
    leg = {'data','shifted'};
    xNull = sort(null, 1, 'ascend');
    xNull = sort(xNull, 2, 'ascend');
    limNull = prctile(xNull, [2.5 97.5], 2);
    limNull = [[1 1].*xLimits(1); limNull; [1 1].*xLimits(2)];
    xNull = median(xNull, 2);
    yNull = (1:length(xNull))' ./ (length(xNull));
    xNull = [xLimits(1); xNull; xLimits(2)];
    yNull = [0; yNull; 1];
    
    h(2) = fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', colors(2,:), ...
        'FaceAlpha', alpha);
end

x = sort(data);
y = (1:length(x))' ./ length(x);
x = [xLimits(1); x; xLimits(2)];
y = [0; y; 1];

h(1) = plot(x, y, 'Color', colors(1,:), 'LineWidth', 2);
legend(h, leg, 'Location', 'NorthWest')