 function [xs, info] = pwls_admm2(x, A, W, yi, C, varargin)
%function [xs, info] = pwls_admm2(x, A, W, yi, C, [options])
%|
%| penalized weighted least squares (PWLS) image reconstruction/restoration
%| with convex non-quadratic regularization,
%| cost(x) = (y-Ax)'W(y-Ax)/2 + reg * pot(C * x)
%| minimized via ADMM algorithm with two splits, ala ramani:12:asb
%|	http://dx.doi.org/10.1109/TMI.2011.2175233
%| u = A x and v = C x
%| for which the AL is AL(x,u,v) = (y-u)'W(y-u)/2 + reg * pot(v) ...
%|	+ mu_a/2 |u - A x - eta_u|_2^2 + mu_c/2 |v - C x - eta_v|_2^2
%|
%| in
%|	x	[(N)]		initial estimate
%|	A	[nd (*N)]	system matrix
%|	W	[nd nd]		data weighting matrix, usually Gdiag(wi) or 1
%|	yi	[(nd) 1]	noisy data
%|	C			regularizer matrix (often from Cdiffs)
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

if nargin == 1 && streq(x, 'test'), pwls_admm2_test, return, end
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
arg.mu_a = [];
arg.mu_c = [];
arg.cond_aa_cc = 20; % target condition number for dftAA + factor * dftCC
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
if isempty(arg.mu_a) || isempty(arg.mu_c) || isempty(arg.precon)
	t = num2cell(floor(siz/2+1)); % "center"
	ee = zeros(siz);
	ee(t{:}) = 1;
	AAe = reshape(A' * (A * ee), siz);
	CCe = reshape(C' * (C * ee), siz);
%	im(AAe)
%	im(CCe)
	dftAA = fftn(ifftshift(AAe));
	dftAA = max(reale(dftAA), 0); % A'A is positive semidefinite
	dftCC = fftn(ifftshift(CCe));
	dftCC = max(reale(dftCC), 0); % C'C is positive semidefinite
%	minmax(dftAA)
%	minmax(dftCC)
end

% select AL penalty parameters

if isempty(arg.mu_c)
	arg.mu_c = arg.reg; % 'todo' for now
end

if isempty(arg.mu_a)
	b_best = pwls_admm2_best_b(dftAA, dftCC, ... % best b = mu_c / mu_a
			'cond_aa_cc', arg.cond_aa_cc, 'chat', 0);
	arg.mu_a = arg.mu_c / b_best;
	arg.mu_a = 1; % todo: works better!?  still needs more work...
	pr '[b_best arg.mu_c / arg.mu_a]'
end

% preconditioner for (mu1 A'A + mu2 C'C)^-1
if isempty(arg.precon)
	precon_filter = arg.mu_a * dftAA + arg.mu_c * dftCC;
	if any(precon_filter(:) <= 0), 'fail bug', end
	printm('cond # = %g', max(abs(precon_filter)) / min(abs(precon_filter)))
	precon_filter = 1 ./ precon_filter;
	mask = true(siz); % todo: need actual mask!
	precon_fun = @(arg,x) ifftn(fftn(embed(x,mask)) .* precon_filter);
	arg.precon = fatrix2('forw', precon_fun, 'idim', siz, 'odim', siz);

	if 0
		tmp = arg.mu_a * AAe + arg.mu_c * CCe;
		tmp = arg.precon * tmp;
		clf, im(tmp)
		keyboard
	end
end

if isscalar(W)
	inv_W_I = 1 / (W + arg.mu_a);
	printm('W cond # = %g', max(inv_W_I) / min(inv_W_I))
else
	tmp = W * ones(size(yi));
	tmp = 1 ./ (tmp + arg.mu_a);
	printm('W cond # = %g', max(tmp) / min(tmp))
	inv_W_I = Gdiag(tmp);
	clear tmp
end

printm 'ignore the next warning about mismatched odims:'
As = [sqrt(arg.mu_a) * A; sqrt(arg.mu_c) * C];

Ax = A * x;
Cx = C * x; % could be large memory for problems such as helical CT
u = Ax;
v = Cx;
eta_u = 0;
eta_v = 0;

% iterate
ticker(mfilename, 1, arg.niter)
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

if iter > 1
	% update x by PCG, for AL cost function:
	% mu_a/2 |u - A x - eta_u|_2^2 + mu_c/2 |v - C x - eta_v|_2^2
	xold = x;
	ytmp = [arg.mu_a * col(u - eta_u); arg.mu_c * col(v - eta_v)];
	x = qpwls_pcg1(xold(:), As, 1, ytmp, 0, 'precon', arg.precon, ...
		'isave', 'last', arg.arg_pcg{:});
	x = reshape(x, size(xold));
end

	Ax = A * x;
	Cx = C * x;

	% update u:
	% (y-u)'W(y-u)/2 + mu_a/2 |u - A x - eta_u|_2^2
	u = inv_W_I * (W * yi + arg.mu_a * (Ax + eta_u));

	% update v:
	% reg * pot(v) + mu_c/2 |v - C x - eta_v|_2^2
	tmp = Cx + eta_v;
	v = arg.pot.shrink(Cx + eta_v, arg.reg / arg.mu_c);
%	im plc 1 3, im(1, tmp), im(2, v), im(3, v-tmp), keyboard

	% update eta:
	eta_u = eta_u - (u - Ax);
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


% pwls_admm2_best_b()
% find best "regularization parameter" for A'A + b C'C, where b = mu_c/mu_a
function [b_best c_best] = pwls_admm2_best_b(dftAA, dftCC, varargin)

arg.cond_aa_cc = 20;
arg.lb_min = -10;
arg.lb_max = 5;
arg.lb_n = 101;
arg.chat = 0;
arg = vararg_pair(arg, varargin);

cond_fun1 = @(b) ... % condition of A'A + b C'C, where b = mu_c/mu_a
	max(col(dftAA + b * dftCC)) / min(col(dftAA + b * dftCC));

if 1
	nb = arg.lb_n;
	blist = logspace(arg.lb_min, arg.lb_max, nb);
	clist = zeros(1, nb);
	for ii = 1:numel(blist)
		clist(ii) = cond_fun1(blist(ii));
	end

	ibest = imin(clist);
	if ibest == 1 || ibest == nb
		clf, loglog(blist, clist)
		fail('todo: end case, so use cond_aa_cc')
%		cond_fun2 = @(b) (cond_fun1(b) - cond_aa_cc)^2;
%		b_best = fminbnd(cond_fun2, 1e-10, 1e10);
	end

%	tic
	b_best = fminbnd(cond_fun1, blist(ibest-1), blist(ibest+1));
	c_best = cond_fun1(b_best);
	pr c_best
%	toc

	[blist ind] = sort([blist b_best]);
	tmp = [clist c_best];
	clist = tmp(ind);
	if arg.chat
		clf, loglog(blist, clist, '.-', b_best, c_best, 'o')
		axis([min(blist) max(blist) 1 10*clist(1)])
		prompt
	end
end



% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default(x, iter, varargin)
out = [cpu('etoc')];


% pwls_admm2_test
% basic test using a tiny image restoration problem
function pwls_admm2_test

mask = true([8 7]);
% mask(1) = false; % todo mask
A = Gblur(mask, 'psf', ones(3)/9);
xtrue = zeros(size(mask), 'single');
xtrue(end/2, round(end/2)) = 1;
clim = [-0.1 1];

y = A * xtrue;
C = Cdiffs(size(mask));
W = 1;
reg = 2^-12;

im plc 2 3
im(1, xtrue, clim)
im(2, y)

if 1 % analytical solution as reference (because quadratic)
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

pot = potential_fun('quad');
xadmm = pwls_admm2(xinit, A, W, y, C, 'reg', reg, 'niter', 90, ...
	'pot', pot, 'isave', 'last', 'chat', 1);

im(5, xadmm, clim)
im(6, xadmm - xhat, 'xadmm-xhat'), cbar
equivs(xadmm, xhat, 'thresh', 1e-5)
