 function [xs, info] = qpwls_pcg1(x, A, W, yi, C, varargin)
%function [xs, info] = qpwls_pcg1(x, A, W, yi, C, [options])
%|
%| quadratic penalized weighted least squares (QPWLS) via
%| preconditioned conjugate gradients (PCG) algorithm
%| cost(x) = (y-Ax)'W(y-Ax) / 2 + x'C'Cx / 2
%|
%| in
%|	x	[np 1]		initial estimate
%|	A	[nd np]		system matrix
%|	W	[nd nd]		data weighting matrix, usually Gdiag(wi)
%|	yi	[nd 1]		noisy data
%|	C	[nc np]		penalty 'differencing matrix' (0 for unregularized)
%|
%| options
%|	niter			# total iterations (default: 1)
%|					(max # if tol used)
%|	isave	[]		list of iterations to archive (default: 'last')
%|	userfun	@		user defined function handle (see default below)
%|					taking arguments (x, userarg{:})
%|	userarg {}		user arguments to userfun (default {})
%|	precon	[np np]		preconditioner (matrix or object) (or 1)
%|	dircheck 0|1		check descent direction? (default: 1)
%|				set to 0 to save time, if you dare...
%|	stop_diff_tol		stop iterations if norm(xnew-xold)/norm(xnew)
%|				is less than this unitless value.  default: 0
%|	stop_diff_norm		use norm(.,type) for stop rule
%|				choices: 1 | 2 (default) | inf
%|	stop_grad_tol		stop if norm(grad) / y'W y < tol; default: 0
%|	stop_grad_norm		which norm(grad) to use.  default: 2
%|	chat	0|1		verbosity (default 0)
%|
%| out
%|	xs	[np niter]	estimates each iteration
%|	info	[niter 3]	gamma, step size, time each iteration
%|
%| Copyright Jan 1998, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), qpwls_pcg1_test0, return, end
if nargin < 5, help(mfilename), error(mfilename), end

% defaults
arg.precon = 1;
arg.niter = 1;
arg.isave = [];
arg.userfun = @userfun_default;
arg.userarg = {};
arg.key = 1;
arg.dircheck = true; % default is to check descent direction
arg.stop_diff_tol = 0;
arg.stop_diff_norm = 2;
arg.stop_grad_tol = 0;
arg.stop_grad_norm = 2;
arg.chat = 0;

arg = vararg_pair(arg, varargin, 'subs', ...
{'stop_threshold', 'stop_diff_tol'; 'stop_norm_type', 'stop_diff_norm'});

arg.isave = iter_saver(arg.isave, arg.niter);
if arg.stop_diff_tol
	norm_diff = @(x) norm(x, arg.stop_diff_norm);
end
if arg.stop_grad_tol
	norm_grad = @(g) norm(g, arg.stop_grad_norm) / reale(yi' * (W * yi));
end

if ~isreal(yi) && ~isequal(arg.precon, 1)
	persistent warned
	if isempty(warned), warned = 0; end
	if ~warned
		warning 'not 100% sure about the complex preconditioned case'
		warned = 1;
	end
end

cpu etic
if isempty(x), x = zeros(ncol(A),1); end

x = x(:);
np = length(x);
xs = zeros(np, length(arg.isave), 'single');
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = single(x);
end

%info = zeros(arg.niter, ?); % trick: do not initialize because size may change

% initialize projections
ticker(mfilename, 1, arg.niter)
Ax = A * x;
Cx = C * x;

% iterate
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	% (negative) gradient
	ngrad = A' * (W * (yi-Ax)) - C' * Cx;

	if arg.stop_grad_tol && norm_grad(ngrad) < arg.stop_grad_tol
		if arg.chat
			printm('stop at iteration %d with grad %g < %g', ...
				iter, norm_grad(ngrad), arg.stop_grad_tol)
		end
		if isequal(arg.isave, arg.niter) % saving last iterate only?
			xs = single(x); % save 'final' iterate
		else % saving many iterates?
			xs(:, arg.isave > iter) = []; % clear out unused
		end
	return
	end

	% preconditioned gradient
	pregrad = arg.precon * ngrad;

	% search direction
	newinprod = ngrad' * pregrad;
	% fix: should i take the real part?
	newinprod = reale(newinprod, 'warn', 'inprod');
	if iter == 1
		ddir = pregrad;
		gamma = 0;
	else
		if oldinprod == 0
			warn 'inprod=0. going nowhere!'
			gamma = 0;
		else
			gamma = newinprod / oldinprod;	% Fletcher-Reeves
%			gamma = (newinprod - oldgrad' * pregrad) / oldinprod;
		end
		ddir = pregrad + gamma * ddir;
	end
	oldgrad = ngrad;
	oldinprod = newinprod;

	% check if descent direction
	if arg.dircheck && real(ddir' * ngrad) < 0
		warn 'wrong direction; try using stop_grad_tol'
		ratio = norm(ngrad(:), arg.stop_grad_norm) / (yi'*W*yi);
		pr ratio % see how small it is
		if arg.key, keyboard, end
	end

	% step size in search direction
	Adir = A * ddir;
	Cdir = C * ddir;

	denom = Adir'*(W*Adir) + Cdir'*Cdir;
	denom = reale(denom, 'error', 'denom');
	if denom == 0
		warning 'found exact solution??? step=0 now!?'
		step = 0;
	else
		step = (ddir' * ngrad) / denom;
%		step = reale(step, 'warn', 'step');
		step = real(step); % real step sizes seems only logical
	end

	if step < 0
		warning 'downhill?'
		if arg.key, keyboard, end
	end

	% update
	Ax = Ax + step * Adir;
	Cx = Cx + step * Cdir;
	x = x + step * ddir;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = single(x);
	end
	info(iter,:) = arg.userfun(x, arg.userarg{:});

	% check norm(xnew-xold) / norm(xnew) vs threshold
	if arg.stop_diff_tol && ...
		norm_diff(step * ddir) / norm_diff(x) < arg.stop_diff_tol
		if arg.chat
			ratio = norm_diff(step * ddir) / norm_diff(x);
			printm('stop at iteration %d with diff %g < %g', ...
				iter, ratio, arg.stop_diff_tol)
		end
		if isequal(arg.isave, arg.niter) % saving last iterate only?
			xs = single(x); % save the 'final' iterate
		else % saving many iterates?
			xs(:, arg.isave > iter) = []; % clear out unused
		end
	return
	end
end


% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default(x, varargin)
gamma = evalin('caller', 'gamma');
step = evalin('caller', 'step');
out = [gamma step cpu('etoc')];


% qpwls_pcg1_test0()
function qpwls_pcg1_test0
mask = true([8 7]); mask(1) = false;
A = Gblur(mask, 'psf', ones(3)/9);
tmp = ones(size(mask));
A = Gdiag(tmp(mask), 'mask', mask);
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
im clf, im(xhat)

% user functions for tracking time and distance to a reference image
f.userfun = @(x, xref) [cpu('etoc') norm(x(:) - xref(:))];
f.userarg = {xhat(mask)}; % reference image just for testing

xinit = 0 * mask;
xpcg = qpwls_pcg1(xinit(mask(:)), A, 1, y, R.C, 'niter', 100, ...
	'userfun', f.userfun, 'userarg', f.userarg, ...
        'stop_grad_tol', 1e-8, 'stop_grad_norm', 2, ...
        'stop_diff_tol', 0e-6, 'stop_diff_norm', 2, 'chat', 1);
xpcg = embed(xpcg, mask);

im plc 2 2
im(1, xtrue)
im(2, xhat)
im(3, xpcg)
im(4, xpcg - xhat)
equivs(xpcg, xhat)
