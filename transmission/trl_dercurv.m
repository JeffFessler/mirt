 function [deriv, curv] = trl_dercurv(data, li, curvtype, iblock, nblock)
%function [deriv, curv] = trl_dercurv(data, li, curvtype, iblock, nblock)
%
% evaluate derivatives and curvatures for monoenergetic Poisson transmission
% negative log-likelihoods: f(x) = \sum_i h_i(l), h_i(l) = mi(l) - yi log mi(l)
% mi(l) = bi exp(-l) + ri
%
% in
%	data	{yi, bi, ri}
%				yi is [nb,na].  bi, ri are same or scalar
%				passed to trl_curvature()
%	li	[nb,#view_in_block]
%	curvtype ''		curvature type
%	iblock, nblock		for OS type methods
%
% out
%	deriv	[nb,#]		\dot hi(l)
%	curv	[nb,#]		surrogate curvature for hi at li
%
% Copyright 2004-2-1, Jeff Fessler, The University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end
if nargin < 3, curvtype = 'pc'; end

yi = data{1};
bi = data{2};
ri = data{3};

if nargin == 5
	ia = iblock:nblock:size(yi,2);
	yi = yi(:,ia);
	if length(bi) > 1
		bi = bi(:,ia);
	end
	if length(ri) > 1
		ri = ri(:,ia);
	end
end


% transmission Poisson likelihood function
ei = exp(-li);
mi = bi .* ei + ri;
deriv = (1 - yi ./ mi) .* (-bi .* ei);

curv = trl_curvature(yi, bi, ri, li, curvtype);
