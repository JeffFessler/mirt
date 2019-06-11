 function [xfbp, sino] = em_fbp(sg, ig, yi, ci, ri, varargin)
%function [xfbp, sino] = em_fbp(sg, ig, yi, ci, ri, [options])
% Emission FBP reconstruction from Poisson measurements
% model: Y_i ~ Poisson(c_i [G x]_i + r_i)
% in
%	sg		sino_geom()
%	ig		image_geom()
%	yi		transmission sinogram
%	ci		calibration factors
%	ri		background (randoms, scatter, crosstalk, etc)
%	ci,ri:		optional (can use empty matrices)
%	yi,ci,ri	must have identical dimensions
% option
%	'kernel'	apodization filter kernel (default: [1/3 1/3 1/3])
% out
%	x [np]		image estimate
%	sino [nb,na]	filtered sinogram
%
% Copyright Apr 2000, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

if ~isvar('ci') || isempty(ci)
	ci = sg.ones;
end
if ~isvar('ri') || isempty(ri)
	ri = sg.zeros;
end

arg.kernel = ones(3,1)/3;
arg = vararg_pair(arg, varargin);

eml_check(yi, ci, ri, 'fbp');
[nb na nz] = size(yi);

tmp = fbp2(sg, ig);
xfbp = ig.zeros('nz', nz);
for iz=1:nz
	proj = (yi(:,:,iz) - ri(:,:,iz)) ./ ci(:,:,iz);
	if any(size(arg.kernel) ~= 1)
		proj = conv2(proj, arg.kernel, 'same'); % filter
	end
	xfbp(:,:,iz) = fbp2(proj, tmp);
	xfbp(:,:,iz) = xfbp(:,:,iz) .* ig.mask;
end
