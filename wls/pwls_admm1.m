  function [xs, info] = pwls_admm1(x, A, W, yi, C, varargin)
%|function [xs, info] = pwls_admm1(x, A, W, yi, C, [options])
%|
%| penalized weighted least squares (PWLS)
%| with convex non-quadratic regularization,
%| cost(x) = (y-Ax)'W(y-Ax)/2 + reg * pot(C * x)
%| minimized via ADMM algorithm with one splits, ala "split Bregman"
%| v = C x
%| for which the AL is AL(x,u,v) = (y-Ax)'W(y-Ax)/2 + reg * pot(v) ...
%|	+ mu_c/2 |v - C x - eta_v|_2^2
%| assumes that we can invert A'WA + C'C with backslash, i.e., sparse!
%|
%| in
%|	x	[(N)]		initial estimate
%|	A	[nd (*N)]	system matrix (sparse)
%|	W	[nd nd]		data weighting matrix, usually sparse or 1
%|	yi	[(nd) 1]	noisy data
%|	C			regularizer matrix (often from Cdiffs) sparse
%|
%| options
%|	reg			regularization parameter
%|				(scalar or vector of size C*x)
%|	pot			potential function, see potential_fun.m
%|					default: potential_fun('l1');
%|	niter			# total iterations (default: 1)
%|					(max # if tol used)
%|	isave	[]		list of iterations to archive (default: 'last')
%|	userfun	@		user defined function handle (see default below)
%|					taking arguments (x, iter, userarg{:})
%|	userarg {}		user arguments to userfun (default {})
%|	precon	[np np]		preconditioner
%|				(default: FFT-based for mu1 A'A + mu2 C'C)
%|	stop_diff_tol		stop iterations if norm(xnew-xold)/norm(xnew)
%|				is less than this unitless value.  default: 1e-6
%|	stop_diff_norm		use norm(.,type) for stop rule
%|				choices: 1 | 2 (default) | inf
%|	chat	0|1		verbosity (default 0)
%|
%| out
%|	xs	[(N) niter]	estimates each iteration
%|	info	[niter ?]	time each iteration
%|
%| Copyright 2014-04-27, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), pwls_admm1_test, return, end
if nargin < 5, help(mfilename), error(mfilename), end

warn('todo: under construction')

% defaults
arg.reg = 1;
arg.pot = [];
arg.precon = [];
arg.niter = 1;
arg.isave = [];
arg.userfun = @userfun_default;
arg.userarg = {};
%arg.key = 1;
arg.stop_diff_tol = 1e-6;
arg.stop_diff_norm = 2;
arg.mu_c = [];
arg.cond_awa_cc = 20; % target condition number for dftAWA + factor * dftCC
arg.arg_pcg = {'niter', 2};
arg.chat = 0;

arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);
if arg.stop_diff_tol
	norm_diff = @(x) norm(x(:), arg.stop_diff_norm);
end

cpu etic
if isempty(x), x = zeros(ncol(A),1); end

np = numel(x);
xs = zeros(np, numel(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x(:);
end

%info = zeros(arg.niter, ?); % trick: do not initialize because size may change
if isempty(arg.pot)
	arg.pot = potential_fun('l1');
end

siz = size(x);
if isempty(arg.mu_c) || isempty(arg.precon)
	t = num2cell(floor(siz/2+1)); % "center"
	ee = zeros(siz);
	ee(t{:}) = 1;
	AWAe = reshape(A' * (W * (A * ee(:))), siz);
	CCe = reshape(C' * (C * ee(:)), siz);
%	im(AWAe)
%	im(CCe)
	dftAWA = fftn(ifftshift(AWAe));
	dftAWA = max(reale(dftAWA, 'warn'), 0); % A'A is positive semidefinite
	dftCC = fftn(ifftshift(CCe));
	dftCC = max(reale(dftCC), 0); % C'C is positive semidefinite
%	minmax(dftAWA)
%	minmax(dftCC)
end

% select AL penalty parameter

if isempty(arg.mu_c)
	c_best = pwls_admm1_best_c(dftAWA, dftCC, ... % best b = mu_c / mu_a
			'cond_awa_cc', arg.cond_awa_cc, 'chat', 0);
	arg.mu_c = arg.reg; % 'todo' for now
%	arg.mu_c = c_best; % 'todo' for now
end

% sparse matrix A'WA + mu_c C'C
Hc = A' * (W * A) + arg.mu_c * (C' * C);
%pr cond(full(Hc))

Cx = C * x(:); % could be large memory for problems such as helical CT
v = Cx;
eta_v = 0;

% iterate
ticker(mfilename, 1, arg.niter)
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

if iter > 1
	% update x by PCG, for AL cost function:
	% 1/2 |y - A x|_W^2 + mu_c/2 |v - C x - eta_v|_2^2
	xold = x;
	x = Hc \ (A' * (W * yi(:)) + arg.mu_c * (C' * (v(:) - eta_v(:))));
	x = reshape(x, size(xold));

	Cx = C * x(:);
end

	% update v:
	% reg * pot(v) + mu_c/2 |v - C x - eta_v|_2^2
	tmp = Cx + eta_v;
	v = arg.pot.shrink(Cx + eta_v, arg.reg / arg.mu_c);
%	im plc 1 3, im(1, tmp), im(2, v), im(3, v-tmp), keyboard

	% update eta:
	eta_v = eta_v - (v - Cx);

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x(:);
	end
	info(iter,:) = arg.userfun(x, iter, arg.userarg{:});

	% check norm(xnew-xold) / norm(xnew) vs threshold

	if iter > 1 && arg.stop_diff_tol && ...
		norm_diff(xold-x) / norm_diff(x) < arg.stop_diff_tol
		if arg.chat
			ratio = norm_diff(xold-x) / norm_diff(x);
			printm('stop at iteration %d with diff %g < %g', ...
				iter, ratio, arg.stop_diff_tol)
		end
		if isequal(arg.isave, arg.niter) % saving last iterate only?
			xs = x; % save the 'final' iterate
		else % saving many iterates?
			xs(:, arg.isave > iter) = []; % clear out unused
		end
	return
	end
end
xs = reshapee(xs, size(x), []);


% pwls_admm1_best_c()
% find best "regularization parameter" for A'A + mu_c C'C
function [c_best k_best] = pwls_admm1_best_c(dftAWA, dftCC, varargin)

arg.cond_awa_cc = 20;
arg.lc_min = -10;
arg.lc_max = 5;
arg.lc_n = 101;
arg.chat = 0;
arg = vararg_pair(arg, varargin);

cond_fun1 = @(c) ... % condition of A'WA + mu_c C'C
	max(col(dftAWA + c * dftCC)) / min(col(dftAWA + c * dftCC));

if 1
	nc = arg.lc_n;
	clist = logspace(arg.lc_min, arg.lc_max, nc);
	klist = zeros(1, nc);
	for ii = 1:numel(clist)
		klist(ii) = cond_fun1(clist(ii));
	end

	ibest = imin(klist);
	if ibest == 1 || ibest == nc
		clf, loglog(clist, klist)
		fail('todo: end case, so use cond_awa_cc')
%		cond_fun2 = @(c) (cond_fun1(c) - cond_awa_cc)^2;
%		c_best = fminbnd(cond_fun2, 1e-10, 1e10);
	end

%	tic
	c_best = fminbnd(cond_fun1, clist(ibest-1), clist(ibest+1));
	k_best = cond_fun1(c_best);
	pr k_best
%	toc

	[clist ind] = sort([clist c_best]);
	tmp = [klist k_best];
	klist = tmp(ind);
	if arg.chat
		clf, loglog(clist, klist, '.-', c_best, k_best, 'o')
		axis([min(clist) max(clist) 1 10*klist(1)])
		prompt
	end
end



% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default(x, iter, varargin)
out = [cpu('etoc')];


% pwls_admm1_test
function pwls_admm1_test

mask = true([8 7]);
% mask(1) = false; % todo mask
A = Gblur(mask, 'psf', ones(3)/9);
xtrue = zeros(size(mask));
xtrue(end/2, round(end/2)) = 1;
clim = [-0.1 1];

y = A * xtrue;
C = Cdiffs(size(mask));
W = 1;
reg = 2^-12;

im plc 3 3
im(1, xtrue, clim)
im(2, y)

if 1
	tmp1 = A' * (W * A);
	tmp2 = reg * (C' * C);
	hess = full(tmp1 + tmp2);
	pr cond(hess)
	xhat = hess \ col(A' * (W * y));
	xhat = embed(xhat, mask);
	im(3, xhat, clim)
	%im(4, full(tmp2))
else
	xhat = nan(size(mask));
end

xinit = ir_wls_init_scale(A, y);
im(4, xinit, clim)

As = sparse(A);
Cs = sparse(C);

pot = potential_fun('quad');
xadmm = pwls_admm1(xinit, As, W, y, Cs, 'reg', reg, 'niter', 90, ...
	'pot', pot, 'isave', 'last', 'chat', 1);

im(5, xadmm, clim)
im(6, xadmm - xhat)
equivs(xadmm, xhat, 'thresh', 1e-5)

pot = potential_fun('l1');
reg = 2^-12;
xadmm = pwls_admm1(xinit, As, W, y, Cs, 'reg', reg, 'niter', 90, ...
	'pot', pot, 'isave', 'last', 'chat', 1);

im(8, xadmm, clim)
xlabelf 'l1'
