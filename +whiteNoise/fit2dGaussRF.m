

function [x, fitResp] = fit2dGaussRF(resp, makeplots)
% resp = your measure of the responsiveness of the neuron at each point in
% space
% x = [amplitude, xCenter, xStdev, yCenter, yStdev, rotation, baseline]
%

    xcoords = 1:size(resp,2);
    ycoords = 1:size(resp,1);

    use2dgWithBsl = 0;

    [x, y] = meshgrid(xcoords, ycoords);
    xdata = zeros(size(x,1), size(x,2), 2);
    xdata(:,:,1) = x;
    xdata(:,:,2) = y;
    
    % x0 is the initial guess
    [maxY, maxX] = find(resp==max(resp(:)),1);
    x0 = [1,maxX,mean(diff(xcoords))*5,maxY,mean(diff(ycoords))*5,0];
    
    if use2dgWithBsl
        
        x0(7) = 0;

        lb = [0,min(xcoords),0,min(ycoords),0,-pi/4, 0];
        ub = [realmax('double'),max(xcoords),(max(xcoords))^2,max(ycoords),(max(ycoords))^2,pi/4, max(resp(:))];
        x = lsqcurvefit(@whiteNoise.D2GaussFunctionRotWithBsl,x0,xdata,resp,lb,ub);

        fitResp = whiteNoise.D2GaussFunctionRotWithBsl(x, xdata);
        
    else
        

        lb = [0,min(xcoords),0,min(ycoords),0,-pi/4];
        ub = [realmax('double'),max(xcoords),(max(xcoords))^2,max(ycoords),(max(ycoords))^2,pi/4];
        x = lsqcurvefit(@whiteNoise.D2GaussFunctionRot,x0,xdata,resp,lb,ub);

        fitResp = whiteNoise.D2GaussFunctionRot(x, xdata);
        
    end
    
    if makeplots
        figure;
        subplot(3,1,1);
        imagesc(xcoords, ycoords, resp);
        colormap hot;
        
        subplot(3,1,2);
        imagesc(xcoords, ycoords, fitResp);
        
        subplot(3,1,3);
        contour(fitResp,1,'Color', 'k', 'LineWidth', 2.0);
        set(gca, 'YDir', 'reverse');
    end
end