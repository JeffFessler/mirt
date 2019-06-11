  function [xs, info] = qpwls_bb1(x, A, W, yi, C, varargin)
%|function [xs, info] = qpwls_bb1(x, A, W, yi, C, [options])
%|
%| quadratic penalized weighted least squares (QPWLS) via
%| preconditioned (todo) Barzilai and Borwein gradient method.
%| cost(x) = (y-Ax)'W(y-Ax) / 2 + x'C'Cx / 2
%|
%| in
%|	x	[np 1]		initial estimate
%|	A	[nd np]		system matrix
%|	W	[nd nd]		data weighting matrix, usually Gdiag(wi) or 1
%|	yi	[nd 1]		noisy data
%|	C	[nc np]		penalty matrix (0 for unregularized)
%|
%| options
%|	niter			# total iterations (default: 100)
%|					(max # if tol used)
%|	isave	[]		list of iterations to archive (default: 'last')
%|	userfun	@		user defined function handle (see default below)
%|					taking arguments (x, userarg{:})
%|	userarg {}		user arguments to userfun (default {})
%|	precon	[np np]		preconditioner (matrix or object) (or 1)
%|	step0			initial step size (default: [], means calculate)
%|	stop_diff_tol		stop iterations if norm(xnew-xold)/norm(xnew)
%|				is less than this unitless value.  default: 0
%|	stop_diff_norm		use norm(.,type) for stop rule
%|				choices: 1 | 2 (default) | inf
%|	stop_grad_tol		stop if norm(grad) / y'W y < tol; default: 1e-4
%|	stop_grad_norm		which norm(grad) to use.  default: 2
%|	chat	0|1		verbosity (default 0)
%|
%| out
%|	xs	[np niter]	estimates each iteration
%|	info	[niter 2]	step size, time each iteration
%|
%| Copyright 2011-07-15, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), qpwls_bb1_test, return, end
if nargin < 5, help(mfilename), error(mfilename), end

% defaults
arg.precon = 1;
arg.niter = 100;
arg.isave = [];
arg.userfun = @userfun_default;
arg.userarg = {};
arg.step0 = [];
arg.stop_diff_tol = 0;
arg.stop_diff_norm = 2;
arg.stop_grad_tol = 1e-4;
arg.stop_grad_norm = 2;
arg.chat = 0;

arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);
if arg.stop_diff_tol
	norm_diff = @(x) norm(x, arg.stop_diff_norm);
end
if arg.stop_grad_tol
	norm_grad = @(g) norm(g, arg.stop_grad_norm) / reale(yi' * (W * yi));
end

if ~isreal(yi)
	warn 'todo: complex case never tested'
end

if ~isequal(arg.precon, 1)
	warn 'todo: precon not tested'
end

if ~isreal(yi) && ~isequal(arg.precon, 1)
	persistent warned
	if isempty(warned), warned = 0; end
	if ~warned
		warning 'todo: not 100% sure about the complex preconditioned case'
		warned = 1;
	end
end

cpu etic
np = size(A,2);
if isempty(x)
	x = zeros(np,1);
else
	x = x(:);
end

xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

%info = zeros(arg.niter, ?); % trick: do not initialize because size may change

% iterate
ticker(mfilename, 1, arg.niter)
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	grad = A' * (W * (A*x - yi)) + C' * C*x; % gradient

	if arg.stop_grad_tol && norm_grad(grad) < arg.stop_grad_tol
		if arg.chat
			printm('stop at iteration %d with grad %g < %g', ...
				iter, norm_grad(grad), arg.stop_grad_tol)
		end
		if isequal(arg.isave, arg.niter) % saving last iterate only?
			xs = x; % save 'final' iterate
		else % saving many iterates?
			xs(:, arg.isave > iter) = []; % clear out unused
		end
	return
	end

	% preconditioned gradient
	pregrad = arg.precon * grad;

	% step size in search direction
	if iter == 1
		if isempty(arg.step0)
			fail 'todo'
		else
			step = arg.step0;
		end
	else
		tmp = x - xold;
		denom = tmp' * (grad - gradold);
%		denom = reale(denom, 'error', 'denom'); % todo?
		if denom == 0
			warn('quit at iteration %d with denom=0!?', iter)
			if isequal(arg.isave, arg.niter) % saving last iterate only?
				xs = x; % save the 'final' iterate
			else % saving many iterates?
				xs(:, arg.isave > iter) = []; % clear out unused
			end
		return
		end
		step = norm(tmp)^2 / denom;
%		step = real(step); % real step sizes seems only logical
	end

	if step < 0
		warning 'downhill; known to occur for BB'
	end

	% update
	xold = x;
	x = x - step * pregrad;
	gradold = pregrad;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
	info(iter,:) = arg.userfun(x, arg.userarg{:});

	% check norm(xnew-xold) / norm(xnew) vs threshold
	if arg.stop_diff_tol && ...
		norm_diff(step * ddir) < arg.stop_diff_tol * norm_diff(x)
		if arg.chat
			ratio = norm_diff(step * ddir) / norm_diff(x);
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


% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default(x, varargin)
step = evalin('caller', 'step');
out = [step cpu('etoc')];


% qpwls_bb1_test()
function qpwls_bb1_test
mask = true([8 7]); mask(1) = false;
A = Gblur(mask, 'psf', ones(3)/9);
%tmp = ones(size(mask));
%A = Gdiag(tmp(mask), 'mask', mask);
xtrue = zeros(size(mask), 'single');
xtrue(end/2, round(end/2)) = 1;
y = A * xtrue(mask);
beta = 2^-7;
beta = 2^-2;
R = Reg1(mask, 'beta', beta, 'order', 1);
qpwls_psf(A, R.C, 1, mask, 1, 'loop', 0);
hess = full(A' * A + R.C' * R.C);
xhat = hess \ (A' * y);
xhat = embed(xhat, mask);
pr fwhm2(xhat)

im plc 2 2
im(1, xtrue), cbar
im(2, xhat), cbar

% user functions for tracking time and distance to a reference image
f.userfun = @(x, xref) [cpu('etoc') norm(x(:) - xref(:))];
f.userarg = {xhat(mask)}; % reference image just for testing

xinit = 0 * mask;
xbb = qpwls_bb1(xinit(mask(:)), A, 1, y, R.C, 'niter', 100, ...
	'step0', 1, ...
	'userfun', f.userfun, 'userarg', f.userarg, ...
        'stop_grad_tol', 1e-8, 'stop_grad_norm', 2, ...
        'stop_diff_tol', 0e-6, 'stop_diff_norm', 2, 'chat', 1);
xbb = embed(xbb, mask);

im(3, xbb), cbar
im(4, xbb - xhat, 'error'), cbar
equivs(xbb, xhat)
