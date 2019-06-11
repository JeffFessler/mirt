 function [deriv, curv] = wls_dercurv(data, li, curvtype, iblock, nblock)
%function [deriv, curv] = wls_dercurv(data, li, curvtype, iblock, nblock)
%
% evaluate derivatives and curvatures for WLS data-fit function:
% f(x) = \sum_i hi(l), hi(l) = 1/2 wi (yi - l)^2
% in
%	data	{yi, wi}	yi is [nb,na].  wi is same, or just "1"
%	li	[nb,#]		where # is the number of "views" in this block
%	curvtype		ignored since irrelevant for WLS
%	iblock, nblock		for OS type methods
% out
%	deriv	[nb,#]		\dot hi(l)
%	curv	[nb,#]		surrogate curvature for hi at li
%
% Copyright 2004-2-1, Jeff Fessler, The University of Michigan

if nargin < 2, help(mfilename), error args, end

yi = data{1};
wi = data{2};
if nargin == 5
	ia = iblock:nblock:size(yi,2);
	yi = yi(:,ia);
	if length(wi) > 1
		wi = wi(:,ia);
	end
end

deriv = wi .* (li - yi);

if length(wi) == 1
	curv = ones(size(deriv));
else
	curv = wi;
end
