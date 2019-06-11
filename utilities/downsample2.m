  function y = downsample2(x, m, varargin)
%|function y = downsample2(x, m, varargin)
%| downsample by averaging by integer factors
%| m can be a scalar (same factor for both dimensions)
%| or a 2-vector
%| in
%|	x	[nx ny]
%| option
%|	'warn'	0|1	warn if non-integer factor. default: 1
%| out
%|	y	[nx/m ny/m]

if nargin == 1 && streq(x, 'test'), downsample2_test, return, end
if nargin < 2, ir_usage, end

arg.warn = 1;
arg = vararg_pair(arg, varargin);

fun = @(x, d) downsample1(x, d, 'warn', arg.warn);

y = fun(x, m(1));
y = fun(y', m(end))';


function downsample2_test

x = reshape(1:24, [4 6]);
y = downsample2(x, 2);
jf_equal(y, [3.5 11.5 19.5; 5.5 13.5 21.5])
if 1 % big test
	x = zeros(2^12);
	cpu etic
	y = downsample2(x, 2);
	cpu etoc 'downsample2 time:'
end
