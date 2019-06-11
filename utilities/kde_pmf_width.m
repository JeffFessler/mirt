 function dx = kde_pmf_width(x, varargin)
%function dx = kde_pmf_width(x, varargin)
%|
%| Determine bin width for kde_pmf1,2
%|
%| in
%|	x	[N 1]	iid realizations of random variable pairs (x,y)
%|
%| option
%|	'type'		'silverman' to choose dx based on rule-of-thumb
%|			using inter-quartile range (p48 of 1986 book)
%|	'chat'	0|1	print out intermediate values? (default: 0)
%|
%| out
%|	dx	[1 1]	bin width
%|
%| Default is for a quadratic spline kernel supported on (-3/2,3/2).
%|
%| Copyright 2010-07-31, Jeff Fessler, University of Michigan

arg.type = 'silverman';
arg.kernel = 'bspline2';
arg.chat = 0;
arg = vararg_pair(arg, varargin);

if nargin < 1, ir_usage, end
if streq(x, 'test'), kde_pmf_width_test, return, end

if ~streq(arg.kernel, 'bspline2')
	fail('only bspline2 done')
end

switch arg.type
case 'silverman'
	dx = kde_pmf_width_silverman(x(:), arg.chat);
otherwise
	fail('not done')
end


% kde_pmf_width_silverman()
function dx = kde_pmf_width_silverman(x, chat)
N = numel(x);
pn = jf_protected_names;
sig = std(x);
iqr = diff(pn.prctile(x, [25 75])); % inter-quartile range

if sig == 0
	warn 'sig = 0 so dx = 0'
	dx = 0;
return
end

if iqr == 0
	warn 'iqr = 0, so ignoring iqr'
	iqr = inf;
end

hx = 0.9 * min(sig, iqr/1.34) / N^(1/5); % rule of thumb for gaussian
dx = hx * sqrt(log(256)) / (3 - sqrt(3)); % convert to bspline2 width

if chat > 1
	pr sig
	pr iqr
end
if chat
	pr dx
end


% test routine
function kde_pmf_width_test

rng(0)
n = 300;
sig = 5;
x = sig * randn(n, 1);
dx = kde_pmf_width(x, 'chat', 2 * im());
