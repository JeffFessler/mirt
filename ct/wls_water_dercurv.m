  function [deriv, curv, fi, di] = wls_water_dercurv(data, li, curvtype, iblock, nblock)
%|function [deriv, curv, fi, di] = wls_water_dercurv(data, li, curvtype, iblock, nblock)
%|
%| evaluate derivatives and curvatures for WLS data-fit function:
%| f(x) = \sum_i hi(l), hi(l) = 1/2 wi (yi - g(l))^2
%| where g() is associated with polyenergetic beam-hardening
%| for a single material type (usually water).
%| Here, g(l) = ftab(scale_x * l)
%|
%| in
%|	data	{yi, wi, ftab, scale_x}
%|				yi is [nb,na].  wi is same, or just "1"
%|	li	[nb,#]		where # is the number of "views" in this block
%|	curvtype		ignored since irrelevant for WLS
%|	iblock, nblock		for OS type methods
%| out
%|	deriv	[nb,#]		\dot hi(l)
%|	curv	[nb,#]		surrogate curvature for hi at li
%|	fi	[nb,#]		g(li), should approximate yi
%|	di	[nb,#]		\dot g(li)
%|
%| todo: this could generalize to any nonlinear functions g!
%|
%| Copyright 2008-9-27, Jeff Fessler

if nargin < 2, help(mfilename), error args, end

[yi wi ftab scale_x] = deal(data{:});

if nargin == 5 % handle subset case
	ia = iblock:nblock:size(yi,2);
	yi = yi(:,ia);
	if length(wi) > 1
		wi = wi(:,ia);
	end
end

fit = ftab.fit;

fi = fit.fmfun(scale_x * li(:));
fi = reshape(fi, size(yi));
di = scale_x * fit.fgrad(scale_x * li(:));
di = reshape(di, size(yi));
ci = scale_x^2 * fit.ls_curv(scale_x * li(:), fi(:));
ci = reshape(ci, size(yi));

deriv = wi .* di .* (fi - yi);
curv = wi .* ci;
