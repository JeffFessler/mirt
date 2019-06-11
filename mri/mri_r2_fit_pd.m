  function pd = mri_r2_fit_pd(images, te, r2, varargin)
%|function pd = mri_r2_fit_pd(images, te, r2, varargin)
%|
%| Estimate pd (proton density) map given R2=1/T2 estimate
%|
%| in
%|	images	[(nd) nt]	images with nt different echo times
%|	te	[nt]		echo times
%|	r2	[(nd)]		estimate of R2 map
%|
%| out
%|	pd	[(nd)]		estimated pd map
%|
%| 2012-05-16, Jeff Fessler, University of Michigan

%if nargin == 1 && streq(x, 'test'), mri_r2_fit_pd_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.chat = 0;
arg = vararg_pair(arg, varargin);

nt = numel(te);
dims = size(images);
if nt ~= dims(end)
	fail('numel(te) = %d but last dimension of images is %d', ...
		numel(te), dims(end))
end

num = 0;
den = 0;
for it=1:nt
	tmp = exp(-te(it) * r2);
	num = num + stackpick(images, it) .* tmp;
	den = den + tmp.^2;
end
pd = num ./ den;
