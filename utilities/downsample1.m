  function y = downsample1(x, down, varargin)
%|function y = downsample1(x, down, varargin)
%| downsample by factor m along first dimension by averaging
%|
%| in
%|	x	[n1 (Nd)]
%|	down			integer downsampling factor
%| option
%|	'warn'	0|1		warn if noninteger multiple (default: 1)
%| out
%|	y	[n1/down (Nd)]
%|
%| Copyright 2009-6-4, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), downsample1_test, return, end
if nargin < 2, ir_usage, end

arg.warn = 1;
arg = vararg_pair(arg, varargin);

dim = size(x);
x = reshape(x, dim(1), prod(dim(2:end))); % [n1 *Nd]
m1 = floor(dim(1) / down);
if m1 * down < dim(1)
	if arg.warn
		warn('truncating input size %d to %d', dim(1), m1 * down)
	end
	x = x(1:(m1*down),:);
end
y = reshapee(x, down, []);
y = mean(y, 1);
y = reshape(y, [m1 dim(2:end)]);


if 0
	% old way with loop
	n1 = floor(size(x,1) / m);
	n2 = size(x,2);
	y = zeros(n1,n2);
	for ii=0:m-1
		y = y + x(ii+[1:m:m*n1],:);
		ticker(mfilename, ii+1,m)
	end
	y = y / m;
end


function downsample1_test

x = reshape(2:2:48, [4 6]);
y = downsample1(x, 2);
jf_equal(y, reshape(3:4:47, [2 6]))
if 0 % big test
	x = zeros(2^12);
	cpu etic
	y = downsample1(x, 2);
	cpu etoc 'downsample1 time:'
end
