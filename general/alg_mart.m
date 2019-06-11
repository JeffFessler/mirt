  function xs = alg_mart(x, At, y, varargin)
%|function xs = alg_mart(x, At, y, [options])
%|
%| multiplicative ART (MART) algorithm.
%| tries to maximize entropy?
%|
%| in
%|	x	[np 1]		initial guess.  must be positive!
%|	At	[np nd]		*transpose* of system matrix
%|	y	[nb na]		measurements
%|
%| option
%|	niter			# of iterations
%|	isave			which iterations to save (default: [] 'last')
%|	eps
%|
%| out
%|	xs	[np niter]	iterates
%|
%| Copyright 2006-4-2, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), alg_art1_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end
if isempty(x), x = ones(nrow(G),1); end

% defaults
arg.niter = 1;
arg.isave = [];
%arg.pixmax = inf;
%arg.pixmin = -inf;

% options 
arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);

[nb na] = size(y);
starts = subset_start(na);

asum = sum(At);

iglist = col(outer_sum(1:nb, (starts-1)*nb));
iglist = iglist(asum(iglist) ~= 0);

np = length(x);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

ticker(mfilename, 1, arg.niter)

for iter=1:arg.niter
	ticker(mfilename, iter, arg.niter)

	for ii=iglist'
%		ai = At(:,ii);
%		x = x .* (y(ii) / (ai' * x)) .^ ai;

%		todo: try following approach to see if faster
		[j ignore ai] = find(At(:,ii));
		xj = x(j);
		inprod = ai' * xj;
		if inprod
			x(j) = xj .* (y(ii) / inprod) .^ ai;
		end
	end

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
end
