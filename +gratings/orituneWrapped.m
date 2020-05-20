function f = orituneWrapped(pars,oris) 
% ORITUNEWRAPPED 	sum of two wrapped gaussians, for orientation tuning
% 
%		orituneWrapped([Dp, Rp, Rn, Ro, sigma], oris), where
%		Dp is the preferred direction (bet 0 and 360)
%		Rp is the response to the preferred direction;
%		Rn is the response to the opposite direction;
%		Ro is the background response (useful only in some cases)
%		sigma is the tuning width;
% 		oris are the orientations, bet 0 and 360.
%
% See also: ORITUNE, FITORI, CIRCSTATS, CIRCSTATS360

% 2015 Sylvia Schroeder (based on oritune)
% part of the Matteobox toolbox

if length(pars(:))~= 5
    error('You must specify all 5 parameters');
end

Dp = pars(1); Rp = pars(2); Rn = pars(3); Ro = pars(4); sigma = pars(5);

siz = size(oris);
oris = reshape(oris,1,[]);
gauss1 = Rp.*sum(exp((-bsxfun(@plus,oris-Dp,(-5:5)'.*360).^2)./...
    (2*sigma^2)));
gauss2 = Rn.*sum(exp((-bsxfun(@plus,oris-Dp+180,(-5:5)'.*360).^2)./...
    (2*sigma^2)));

f = Ro + gauss1 + gauss2;
f = reshape(f, siz);