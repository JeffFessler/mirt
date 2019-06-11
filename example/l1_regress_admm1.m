  function [xs, info] = l1_regress_admm1(yi, A, varargin)
%|function [xs, info] = l1_regress_admm1(yi, A, varargin)
%|
%| l1 regression, minimizing cost function
%| cost(x) = pot( y - A x ) where pot(r) = |r|_1 by default.
%| minimized via ADMM algorithm with v = A x split.
%|
%| in
%|	yi	[M 1]		noisy data
%|	A	[M N]		matrix
%|
%| options
%|	x0	[N 1]		initial estimate (default: A \ yi)
%|	shrink			shrink function for pot: shrink(x, reg)
%|				(default: l1, soft thresholding)
%|
%|	niter			# total iterations (default: 1)
%|					(max # if tol used)
%|	isave	[]		list of iterations to archive (default: 'last')
%|	userfun	@		user defined function handle (see default below)
%|					taking arguments (x, iter, userarg{:})
%|	userarg {}		user arguments to userfun (default {})
%|	stop_diff_tol		stop iterations if norm(xnew-xold)/norm(xnew)
%|				is less than this unitless value.  default: 0
%|	stop_diff_norm		use norm(., type) for stop rule
%|				choices: 1 | 2 (default) | inf
%|	chat	0|1		verbosity (default 0)
%|
%| out
%|	xs	[N niter]	estimates each iteration
%|	info	[niter 1]	time each iteration (for default userfun)
%|
%| Copyright 2013-03-21, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(yi, 'test'), ir_regress_admm1_test, return, end

% defaults
arg.x0 = [];
arg.C = [];
arg.beta = 1; % regularization parameter
arg.shrink = [];
arg.rho = 1; % AL penalty parameter

arg.niter = 1;
arg.isave = [];
arg.userfun = @userfun_default;
arg.userarg = {};
arg.stop_diff_tol = 0;
arg.stop_diff_norm = 2;
arg.chat = 0;

arg = vararg_pair(arg, varargin);

if isempty(arg.x0)
	arg.x0 = A \ yi;
end

x = arg.x0;

if isempty(arg.shrink)
	tmp = potential_fun('l1', 1);
	shrink = @(z, reg) tmp.shrink(z, reg);
else
	shrink = arg.shrink;
end

arg.isave = iter_saver(arg.isave, arg.niter);
if arg.stop_diff_tol
	norm_diff = @(x) norm(x(:), arg.stop_diff_norm);
end

cpu etic

np = numel(x);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x(:);
end

eta = 0; % dual variable

rho = arg.rho;

% iterate
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	xp = x; % previous

	% update primal
	x = A \ (yi + eta); 

	% update auxiliary
	tmp = A * x - eta - yi;
	vv = yi + shrink(tmp, 1 / rho);

	% update multiplier
	tmp = A * x - vv;
	eta = eta - tmp;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
	info(iter,:) = arg.userfun(x, iter, arg.userarg{:});

	% check norm(xnew-xold) / norm(xnew) vs threshold
	if iter > 1 && arg.stop_diff_tol
		ratio = norm_diff(x - xp) / norm_diff(x);
		if ratio < arg.stop_diff_tol
			if 1 || arg.chat
				printm('stop at iteration %d, diff %g < %g', ...
					iter, ratio, arg.stop_diff_tol)
			end
			if isequal(arg.isave, arg.niter) % saving only last?
				xs = x; % save the 'final' iterate
			else % saving many iterates?
				xs(:, arg.isave > iter) = []; % clear out unused
			end
		break
		end
	end
end


% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default(x, iter, varargin)
%pr minmax(x)
out = [cpu('etoc')];
