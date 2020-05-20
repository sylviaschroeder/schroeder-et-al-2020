function [pars, err] = fitoriConditions( oo, rr, cc, paramDeltas, ...
    paramLimits, n, parsInit)
% FITORIWRAPPED fit orientation tuning data with orituneWrapped curve,
% which uses wrapped Gaussians (peak repeats every 360 degrees) instead of 
% simple Gaussians; two curves to two datasets are fit while keeping some
% parameters the same
% 
%		syntax is pars = fitoriTwoConditions(oo,rr), where oo are the orientations,
%		rr are the responses, and pars = [Dp, Rp, Rn, Ro, sigma];
%
%		fitoriWrapped(xx) assumes oo = xx(1,:) and rr = xx(2,:)
% 
%		fitori(oo,rr,ss) lets you specify the weights ss
%		(default: [])
%
%		fitori(oo,rr,ss,fixedpars) lets you specify the value of certain
%		parameters (put NaN if you don't want to specify). DEFAULT: [],
%		which means all NaNs.
%
%		fitoriWrapped(oo,rr,ss,fixedpars,oriflag) if oriflag is set to 'ori' Rn is
%		set to zero (useful if data are in the range of 0 to 180).
%
% Example:
% % underlying parameters (see oritune for their meaning)
% realpars = [30 100 10 5 20]
% 
% % make the data based on those parameters
% oo = 0:15:360;
% dd = orituneWrapped(realpars,oo) + normrnd( 0, 5, 1, length(oo));
% 
% % fit the parameters
% pars = fitoriWrapped(oo,dd);
% 
% % look at the results
% figure; 
% plot(oo,dd, 'ko', 'MarkerFaceC','k'); hold on
% plot( 0:360, orituneWrapped(pars, 0:360), 'k-' );
% plot( 0:360, orituneWrapped(realpars, 0:360), ':' );
% set(gca,'xtick',0:45:360);
%
%       SEE ALSO: orituneWrapped, oritune, circstats, circstats360

% 1997 Matteo Carandini
% 1999 Matteo Carandini, corrected pref ori bias
% 2008-05 LB added oriflag
% 2008-11 LB changed maxRp to the 2*max(rr) instead of max(rr)
% 2009-03 LB fixed a small bug when estimating the start value of sigma
% 2009-03 LB fixed some bugs related to oriflag
% 2009-03 LB added the output argument err
% 2009-05 SK added n
% 2009-07 MC added example of usage
% 2014-06 MC changed behavior of 'ori' flag, used to impose Rp=Rn, but
% that's buggy
% 2015-10 SS introduced wrapped Gaussians (also orituneWrapped)

% part of the Matteobox toolbox

if nargin < 7
    parsInit = [];
end
if nargin < 6
    n = 5;
end
if nargin <5
    paramLimits = repmat([-Inf;Inf],1,5);
end
if nargin < 4 || isempty(paramDeltas)
    paramDeltas = NaN(1,5);
end

oo = oo(:);
rr = rr(:);
cc = cc(:);

if any(size(oo)~=size(rr)) || any(size(oo)~=size(cc))
	error('oo, rr and cc have different sizes');
end

numCond = length(unique(cc));
if ~isempty(parsInit) && length(parsInit) ~= numCond*5
    if length(parsInit) >= 5
        parsInit(6:end) = [];
    else
        parsInit = [];
    end
end

%---------------- make all the rr be>0

minrr = min(rr);
rr = rr-minrr;

%---------------- fit all data to set initial conditions

fixedPars = NaN(1,5);
if isempty(parsInit)
    parsInit = gratings.fitoriWrapped(oo, rr, [], fixedPars, '', n);
    parsInit(3) = (parsInit(2)-parsInit(3))/(parsInit(2)+parsInit(3)); % DI
    parsInit(1) = mod(parsInit(1), 360);
else
    parsInit(4) = parsInit(4) - minrr;
    paramLimits(:,4) = paramLimits(:,4) - minrr;
end

for k = 1:10

    %---------------- preferred ori
    
    prefori = mod(parsInit(1),360);
    minprefori = max(0, paramLimits(1,1));
    maxprefori = min(360, paramLimits(2,1));
    prefori = min(max(prefori,minprefori),maxprefori);
    
    %---------------- resp to pref ori
    
    rp = parsInit(2);
    minrp = max(0, paramLimits(1,2));
    maxrp = min(2*max(rr), paramLimits(2,2));
    rp = min(max(rp,minrp),maxrp);
    
    %---------------- direction index
    
    di = parsInit(3);
    mindi = max(-1, paramLimits(1,3));
    maxdi = min(1, paramLimits(2,3));
    di = min(max(di,mindi),maxdi);
    
    %---------------- offset
    
    r0 = parsInit(4);
    minr0 = max(-max(rr)/5, paramLimits(1,4));
    maxr0 = min(max(rr), paramLimits(2,4));
    r0 = min(max(r0,minr0),maxr0);
    
    %---------------- tuning width sigma
    
    sigma = parsInit(5);
    minsigma = max(0.1*median(diff(unique(oo))), paramLimits(1,5));
    maxsigma = min(180, paramLimits(2,5));
    sigma = min(max(sigma,minsigma),maxsigma);
    
    %---------------- delta preferred ori
    
    if isnan(paramDeltas(1))
        if length(parsInit) > 5
            deltaprefori = parsInit(6:5:end)';
        else
            deltaprefori = zeros(numCond-1, 1);
        end
        mindeltaprefori = max(0, paramLimits(1,1))-prefori;
        maxdeltaprefori = min(360, paramLimits(2,1))-prefori;
    else
        deltaprefori = ones(numCond-1,1) .* paramDeltas(1);
        mindeltaprefori = paramDeltas(1);
        maxdeltaprefori = paramDeltas(1);
    end
    deltaprefori = min(max(deltaprefori,mindeltaprefori),maxdeltaprefori);
    
    %---------------- delta resp to pref ori
    
    if isnan(paramDeltas(2))
        if length(parsInit) > 5
            deltarp = parsInit(7:5:end)';
        else
            deltarp = zeros(numCond-1, 1);
        end
        mindeltarp = max(0, paramLimits(1,2))-rp;
        maxdeltarp = min(2*max(rr), paramLimits(2,2))-rp;
    else
        deltarp = ones(numCond-1,1) .* paramDeltas(2);
        mindeltarp = paramDeltas(2);
        maxdeltarp = paramDeltas(2);
    end
    deltarp = min(max(deltarp,mindeltarp),maxdeltarp);
    
    %---------------- delta direction index
    
    if isnan(paramDeltas(3))
        if length(parsInit) > 5
            deltadi = parsInit(8:5:end)';
        else
            deltadi = zeros(numCond-1, 1);
        end
        mindeltadi = max(-1, paramLimits(1,3))-di;
        maxdeltadi = min(1, paramLimits(2,3))-di;
    else
        deltadi = ones(numCond-1,1) .* paramDeltas(3);
        mindeltadi = paramDeltas(3);
        maxdeltadi = paramDeltas(3);
    end
    deltadi = min(max(deltadi,mindeltadi),maxdeltadi);
    
    %---------------- delta offset
    
    if isnan(paramDeltas(4))
        if length(parsInit) > 5
            deltar0 = parsInit(9:5:end)';
        else
            deltar0 = zeros(numCond-1, 1);
        end
        mindeltar0 = max(-max(rr)/5, paramLimits(1,4))-r0;
        maxdeltar0 = min(max(rr), paramLimits(2,4))-r0;
    else
        deltar0 = ones(numCond-1,1) .* paramDeltas(4);
        mindeltar0 = paramDeltas(4);
        maxdeltar0 = paramDeltas(4);
    end
    deltar0 = min(max(deltar0,mindeltar0),maxdeltar0);
    
    %---------------- delta tuning width sigma
    
    if isnan(paramDeltas(5))
        if length(parsInit) > 5
            deltasigma = parsInit(10:5:end)';
        else
            deltasigma = zeros(numCond-1, 1);
        end
        mindeltasigma = max(0.1*median(diff(unique(oo))), paramLimits(1,5)) - sigma;
        maxdeltasigma = min(180, paramLimits(2,5)) - sigma;
    else
        deltasigma = ones(numCond-1,1) .* paramDeltas(5);
        mindeltasigma = paramDeltas(5);
        maxdeltasigma = paramDeltas(5);
    end
    deltasigma = min(max(deltasigma,mindeltasigma),maxdeltasigma);
    
    %---------------------------- finally, do the fit -------------------------
    
    [ err, pars ] = models.fitit(  'gratings.orituneWrappedConditions', rr,...
        [ minprefori minrp mindi minr0 minsigma ...
          repmat([mindeltaprefori mindeltarp mindeltadi ...
          mindeltar0 mindeltasigma], 1, numCond-1)],...
        [    prefori    rp    di    r0    sigma ...
          [deltaprefori deltarp deltadi deltar0 deltasigma] ],...
        [ maxprefori maxrp maxdi maxr0 maxsigma ...
        repmat([maxdeltaprefori maxdeltarp maxdeltadi ...
        maxdeltar0 maxdeltasigma], 1, numCond-1)],...
        [0 1e-4 1e-4 n 2000 0], oo, cc );
    
    % convert parameters (i.e. R_p and DI) so that DIs for all conditions
    % are positive
    pp = num2cell(pars);
    [Dp, Rp, DI] = deal(pp{1:3});
    if DI < 0
        Dp = mod(Dp+180,360);
        DI = -DI;
        Rn = Rp;
        Rp = (1+DI)/(1-DI)*Rn;
        for c = 1:numCond-1
            pars(1+5*c) = pars(1+5*c) - 180; % deltaDp
            pars(3+5*c) = pars(3+5*c) - 2*DI; % deltaDI
            pars(2+5*c) = Rn + pars(2+5*c) - Rp; % deltaRp
        end
    end
    for c = 1:numCond-1
        if DI+pars(3+5*c) < 0
            pars(1+5*c) = pars(1+5*c) + 180; % deltaDp
            pars(3+5*c) = -pars(3+5*c) - 2*DI; % deltaDI
            Rn2 = Rp+pars(2+5*c);
            Rp2 = (1+DI+pars(3+5*c))/(1-DI-pars(3+5*c)) * Rn2;
            pars(2+5*c) = Rp2 - Rp; % deltaRp
        end
    end
    pars(1:3) = [Dp Rp DI];
    
    % now check whether any parameter limits are unsatisfied; if so, redo the
    % fitting with the current fit as starting point
    if all([pars(2),pars(2) + pars(7:5:end)] >= minrp) && ...
            all([pars(2),pars(2) + pars(7:5:end)] <= maxrp) && ...
            all([pars(3),pars(3)+pars(8:5:end)] <= maxdi) && ...
            all(pars(4) + pars(9:5:end) >= minr0) && ...
            all(pars(4) + pars(9:5:end) <= maxr0) && ...
            all(pars(5) + pars(10:5:end) >= minsigma) && ...
            all(pars(5) + pars(10:5:end) <= maxsigma)
        break
    end
    parsInit = pars;
end
               
% add minrr back to the pars:
pars(4) = pars(4)+minrr;