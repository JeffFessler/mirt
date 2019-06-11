  function [lo hi] = ir_dwt_filters(wname, varargin)
%|function [coef codes] = ir_odwt1(x, varargin)
%|
%| filters for discrete wavelet transform DWT
%|
%| in
%|	'wname'	char	'haar' (default) or 'db4' or 'sym2'
%|
%| option
%|	'usemat' 0|1	default: 0. (if 1 then use matlab wfilters)
%|	'ortho'	0|1	default: 1 (make normalized)
%|	'abs'	0|1	default: 0. (if 1 then take abs of filters)
%|
%| out
%|	lo	[K]	low-pass decomposition filter
%|	hi	[K]	hi-pass decomposition filter
%|
%| 2012-05-21, Jeff Fessler, Univ. of Michigan
%| 2014-01-23, Matt Muckley: added 'db4'

if nargin < 1, ir_usage, end
if streq(wname, 'test'), ir_dwt_filters_test, return, end

arg.wname = wname;
arg.ortho = true;
arg.usemat = false;
arg.abs = false;
arg = vararg_pair(arg, varargin);

if isempty(arg.wname)
	 arg.wname = 'haar';
end

% decomposition filters
if arg.usemat
	[lo hi lo_r hi_r] = wfilters(arg.wname);
	lo = lo'; lo_r = lo_r';
	hi = hi'; hi_r = hi_r';

	jf_equal(lo_r, flipud(lo))
	jf_equal(hi_r, flipud(hi))
else
	switch arg.wname
	case 'haar'
		lo = [1 1]';
		hi = [-1 1]';

	case 'db4'
		h1 = (1 + sqrt(3))/4;
		h2 = (3 + sqrt(3))/4;
		h3 = (3 - sqrt(3))/4;
		h4 = (1 - sqrt(3))/4;
		lo = [h1 h2 h3 h4]';
		hi = [h4 -h3 h2 -h1]';
		clear h1 h2 h3 h4;

	case 'sym2'
		lo = [-0.129409522550921; 0.224143868041857; ...
			0.836516303737469; 0.482962913144690];
		hi = [-0.482962913144690; 0.836516303737469; ...
			-0.224143868041857; -0.129409522550921];
	otherwise
		fail('unknown wname "%s"', arg.wname)
	end
end

if arg.ortho
	lo = lo / norm(lo);
	hi = hi / norm(hi);
end

if arg.abs
	lo = abs(lo);
	hi = abs(hi);
end


% ir_dwt_filters_test()
function ir_dwt_filters_test

list = {'haar', 'sym2', 'db4'};
for ii=1:numel(list)
	wname = list{ii};

	[lo0 hi0] = ir_dwt_filters(wname, 'usemat', 0);
	if numel(lo0) == 2
		jf_equal(lo0' * hi0, 0)
	else
		equivs(1, 1 + lo0' * hi0)
	end

	if exist('dwt', 'file') == 2 && exist('wfilters', 'file') == 2
		[lo1 hi1] = ir_dwt_filters(wname, 'usemat', 1);
	end

	[lo hi] = ir_dwt_filters(wname, 'ortho', 0);
	jf_equal(lo/norm(lo), lo0)
	jf_equal(hi/norm(hi), hi0)

	[lo hi] = ir_dwt_filters(wname, 'abs', 1);
	jf_equal(abs(lo0), lo)
	jf_equal(abs(hi0), hi)
end
