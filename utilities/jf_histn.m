  function [hist center] = jf_histn(data, varargin)
%|function [hist center] = jf_histn(data, varargin)
%|
%| Fast histogram of multidimensional data for equally-spaced bins.
%| todo: use accumarray?
%|
%| in
%|	data	[N M]	data values to be binned (M-dimensional)
%|
%| option
%|	'min'	[M]	minimum bin values for each dimension (left side)
%|	'max'	[M]	maximum bin values for each dimension (right side)
%|	'nbin'	[M]	# of bins for each dimension (default: 100)
%|
%| out
%|	hist	[[ncent]]	histogram values: sum(hist(:)) = N
%|	center	{ncent}		cell array of bin centers for each dimension
%|
%| Copyright 2010-07-31, Jeff Fessler, University of Michigan

if nargin == 1 && streq(data, 'test'), jf_histn_test, return, end
if nargin < 1, ir_usage, end

arg.min = [];
arg.max = [];
arg.nbin = [];
arg.chat = 0;
arg = vararg_pair(arg, varargin);

M = size(data,2);

if isempty(arg.nbin)
	arg.nbin = 100;
end
if numel(arg.nbin) == 1
	arg.nbin = arg.nbin * ones(M,1);
end

for id=1:M
	tmp = data(:,id);

	if isempty(arg.min)
		xmin = min(tmp); 
	else
		xmin = arg.min(id); 
	end

	if isempty(arg.max)
		xmax = max(tmp); 
	else
		xmax = arg.max(id); 
	end

	if xmin == xmax
		if xmin == 0
			xmin = -0.5;
			xmax = +0.5;
		else
			xmin = 0.5 * xmin;
			xmax = 1.5 * xmin;
		end
	end

	K = arg.nbin(id);
	tmp = (tmp - xmin) / (xmax - xmin); % [0,1]
	tmp(tmp < 0) = 0;
	tmp(tmp > 1) = 1;
	tmp = 1 + tmp * (K-1); % [1 K]
	data(:,id) = round(tmp);

	if K == 1
		center{id} = (xmin + xmax) / 2;
	else
		center{id} = linspace(xmin, xmax, K);
	end
end

[hist hcent] = hist_bin_int(data);
for id=1:M
	tmp = hcent{id};
	K = arg.nbin(id);
	if min(tmp) ~= 1 || max(tmp) ~= K
		minmax(tmp)
		fail 'todo'
	end
end

% test routine
function jf_histn_test

rng(0)
n = 1000;
sig = 5;
rho = 0.6; % correlated gaussian
Cov = sig * [1 rho; rho 1];
tmp = sqrtm(Cov)
data = randn(n, 2) * sqrtm(Cov);

nbin = [30 30]
[hist cent] = jf_histn(data, 'nbin', nbin);
[xs ys] = deal(cent{:});
if im
	im plc 1 2
	im(1, xs, ys, hist)
	axis equal
	im subplot 2
	plot(data(:,1), data(:,2), '.')
	axis equal
end
