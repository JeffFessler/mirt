  function [images kappa_typical] = mri_r2_fit_normalize(images, te, varargin)
%|function [images kappa_typical] = mri_r2_fit_normalize(images, te, varargin)
%|
%| Normalize images for R2=1/T2 mapping so that typical "kappa" is unity.
%|
%| in
%|	images	[(nd) nt]	images with nt different echo times
%|	te	[nt]		echo times
%|
%| option
%|	'chat'		0|1	if 1, show final kappa map.  default: 0
%|	'threshold'	[0,1]	use median of values above this fraction of max
%|				default: 0.5
%| out
%|	images	[(nd) nt]	normalized images
%|	kappa_typical		scaling factor
%|
%| 2012-05-16, Jeff Fessler, University of Michigan

%if nargin == 1 && streq(x, 'test'), mri_r2_fit_normalize_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.threshold = 0.5;
arg.chat = 0;
arg = vararg_pair(arg, varargin);

kappa = sqrt(mri_r2_fit_kappa2(images, te));

is_big = kappa > arg.threshold * max(kappa(:));
kappa_typical = median(kappa(is_big));

images = images / kappa_typical;

if arg.chat
	kappa = sqrt(mri_r2_fit_kappa2(images, te));
	im(kappa), cbar
end


% mri_r2_fit_kappa2()
function kappa2 = mri_r2_fit_kappa2(images, te)

nt = numel(te);
dims = size(images);
if nt ~= dims(end)
	fail('numel(te) = %d but last dimension of images is %d', ...
		numel(te), dims(end))
end

kappa2 = 0;
for it=1:nt
	kappa2 = kappa2 + te(it) * abs(stackpick(images, it)).^2;
end
