 function [xfbp, proj] = tr_fbp(sg, ig, yi, bi, ri, varargin)
%function [xfbp, proj] = tr_fbp(sg, ig, yi, bi, ri, [options])
% Transmission FBP reconstruction from Poisson measurements
% model: Y_i ~ Poisson(b_i exp(-[G x]_i) + r_i)
% in:
%	sg	sino_geom()
%	ig	image_geom()
%	yi	transmission sinogram
%	bi	blank scan factors
%	ri	background (randoms, scatter, crosstalk, etc)
%	bi,ri:	optional (can use empty matrices)
%	yi,bi,ri must have identical dimensions
% out:
%	x [np]	image estimate
%
% Copyright Mar 2000, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

arg.kernel = ones(3,1)/3;
arg = vararg_pair(arg, varargin);

if ~isvar('bi') || isempty(bi)
	bi = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end

trl_check(yi, bi, ri);

proj = (yi - ri) ./ bi;
if any(size(arg.kernel) ~= 1)
	proj = proj - 1;	% so that zero padding is reasonable
	proj = conv2(proj, arg.kernel, 'same');	% filter (with zero pad)
	proj = proj + 1;
end
if any(proj(:) <= 0), warning 'negative filtered projections', end
proj(proj <= 0) = 1;
proj = -log(proj); % kludge

tmp = fbp2(sg, ig);
xfbp = fbp2(proj, tmp);
xfbp = xfbp .* ig.mask;
