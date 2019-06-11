  function xs = alg_art1(x, At, y, varargin)
%|function xs = alg_art1(x, At, y, [options])
%| classical ART algorithm, aka, Kaczmarz algorithm; tries to solve y=Ax
%|
%| in
%|	x	[np 1]		initial guess, possibly empty
%|	Caution: x must be in the range of A' for convergence!
%|	At	[np nd]		(hermitian) *transpose* of system matrix
%|	y	[nb na]		measurement
%|
%| option
%|	wi	[nb na]		weights
%|	niter			# of iterations
%|	isave			default: [] 'last'
%|	pixmin
%|	pixmax
%|	anorms	[nd 1]		see below
%|	eps
%|
%| out
%|	xs	[np niter]	iterates
%|
%| Copyright 2006-4-2, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), alg_art1_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end
if isempty(x), x = zeros(nrow(A),1); end

% defaults
arg.niter = 1;
arg.isave = [];
arg.anorms = [];
arg.eps = eps;
%arg.pixmax = inf;
%arg.pixmin = -inf;

% options 
arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);

[nb na] = size(y);
starts = subset_start(na);

% For WLS, premultiply y and postmultiply At by W^{1/2}
%Wh = spdiag(sqrt(wi(:)), 'nowarn');
%y = Wh * y(:);
%At = At * Wh;

% weighted row norms
if isempty(arg.anorms)
	arg.anorms = sum(At.^2); % | e_i' A |^2
end

iglist = col(outer_sum(1:nb, (starts-1)*nb));
iglist = iglist(arg.anorms(iglist) ~= 0);

adenom = arg.anorms + eps; % trick:

np = length(x);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

ticker(mfilename, 1, arg.niter)

for iter=1:arg.niter
	ticker(mfilename, iter, arg.niter)

	for ii=iglist'
		ai = At(:,ii);
		step = (y(ii) - ai' * x) / adenom(ii);
		x = x + ai * step;

%		todo: try following approach to see if faster
%		[j ignore ai] = find(At(:,ii));
%		xj = x(j);
%		step = (y(ii) - ai' * xj) / adenom(ii);
%		x(j) = xj - step * ai;
	end

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
end
