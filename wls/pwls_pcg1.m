 function [xs, info] = pwls_pcg1(x, A, W, yi, R, varargin)
%function [xs, info] = pwls_pcg1(x, A, W, yi, R, [options])
%|
%| penalized weighted least squares (PWLS)
%| with convex non-quadratic regularization,
%| minimized via preconditioned conjugate gradient algorithm.
%| cost(x) = (y-Ax)'W(y-Ax)/2 + R(x)
%|
%| in
%|	x	[np 1]		initial estimate
%|	A	[nd np]		system matrix
%|	W	[nd nd]		data weighting matrix, usually Gdiag(wi)
%|	yi	[nd 1]		noisy data
%|	R			penalty object (see Reg1.m)
%|
%| options
%|	niter			# total iterations (default: 1)
%|					(max # if tol used)
%|	isave	[]		list of iterations to archive (default: 'last')
%|	userfun	@		user defined function handle (see default below)
%|					taking arguments (x, iter, userarg{:})
%|	userarg {}		user arguments to userfun (default {})
%|	precon	[np np]		preconditioner (default: 1, i.e., none)
%|	stepper			method for step-size line search
%|				default: {'qs', 3}
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
%| Copyright 1996-7, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), pwls_pcg1_test, return, end
if nargin < 5, ir_usage, end

% defaults
arg.precon = 1;
arg.niter = 1;
arg.isave = [];
arg.stepper = {'qs', 3}; % quad surr with this # of subiterations
arg.userfun = @userfun_default;
arg.userarg = {};
arg.key = 1;
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
	% todo: the "correct" way is sum(abs(sqrtm(W) * y).^2, 'double')
	norm_grad = @(g) norm(g, arg.stop_grad_norm) / reale(yi' * (W * yi));
end

cpu etic
if isempty(x), x = zeros(ncol(A),1, class(yi)); end

x = x(:);
np = numel(x);
xs = zeros(np, numel(arg.isave), class(x));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

%info = zeros(arg.niter, ?); % trick: do not initialize because size may change

% initialize projections
ticker(mfilename, 1, arg.niter)
Ax = A * x;

oldinprod = 0;

% iterate
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	% (negative) gradient
	ngrad = A' * (W * (yi-Ax));
	pgrad = R.cgrad(R, x);
	ngrad = ngrad - pgrad;

	if arg.stop_grad_tol && norm_grad(ngrad) < arg.stop_grad_tol
		if arg.chat
			printm('stop at iteration %d with grad %g < %g', ...
				iter, norm_grad(ngrad), arg.stop_grad_tol)
		end
		if isequal(arg.isave, arg.niter) % saving last iterate only?
			xs = x; % save 'final' iterate
		else % saving many iterates?
			xs(:, arg.isave > iter) = []; % clear out unused
		end
	return
	end

	% preconditioned gradient
	pregrad = arg.precon * ngrad;

	% search direction
	newinprod = dot_double(conj(ngrad), pregrad);
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
	if real(dot_double(conj(ddir), ngrad)) < 0
		warn('wrong direction at iter=%d; try using stop_grad_tol?', iter)
		ratio = norm(ngrad(:), arg.stop_grad_norm) / (yi'*W*yi);
		pr ratio % see how small it is
		if arg.key, keyboard, end
	end

	% step size in search direction
	Adir = A * ddir;
%	Cdir = R.C * ddir; % this is too big for 3D CT problems

	% one step based on quadratic surrogate for penalty
	if streq(arg.stepper{1}, 'qs1')
%		pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir); % avoid Cdir
%		pdenom = (abs(ddir).^2)' * R.denom(R, x);
		pdenom = dot_double(abs(ddir).^2, R.denom(R, x)); % 2012-07-24
%		ldenom = Adir'*(W*Adir);
		% todo: the "correct" way is sum(abs(sqrtm(W) * Adir).^2, 'double')
		ldenom = dot_double(conj(Adir), W*Adir); % 2012-07-24
		ldenom = reale(ldenom);
		denom = ldenom + pdenom;
		if denom == 0
			warning 'found exact solution???  step=0 now!?'
			step = 0;
		else
			step = real((ddir' * ngrad) / denom);
		end

	% Iteratively minimize the 1D line-search function over step:
	%	1/2 || y - A (x + step*ddir) ||_W^2 + R(x + step*ddir)
	% This is a real-valued function of the real-valued step parameter.
	elseif streq(arg.stepper{1}, 'qs')
		nsub = arg.stepper{2};
	%	dAWAd = Adir' * (W * Adir);
		% todo: the "correct" way is sum(abs(sqrtm(W) * Adir).^2, 'double')
		dAWAd = dot_double(conj(Adir), W * Adir); % 2012-07-24
		dAWAd = reale(dAWAd); % 2008-10-16
		dAWr = Adir' * (W * (yi-Ax));
		dAWr = real(dAWr); % 2008-10-16
		step = 0;
		for is=1:nsub
%			pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir); % avoid Cdir
%			pdenom = (abs(ddir).^2)' * R.denom(R, x + step * ddir);
			pdenom = dot_double(abs(ddir).^2, R.denom(R, x + step * ddir));
			denom = dAWAd + pdenom;
			if denom == 0 || isinf(denom)
				if norm(pregrad) == 0
					printm 'found exact solution!?'
					denom = inf; % lazy trick so step=0
				else
					printm '0 or inf denom?'
					if arg.key, keyboard, end
					error bad
				end
			end
			pgrad = R.cgrad(R, x + step * ddir);
			pdot = dot_double(conj(ddir), pgrad);
			pdot = real(pdot); % 2008-10-15
			step = step - (-dAWr + step * dAWAd + pdot) / denom;
%			2008-10-16: removed below because made real above
%			step = real(step); % real step size seems logical
		end

	else
		error 'bad stepper'
	end

	if step < 0
		warning 'downhill?'
		if arg.key, keyboard, end
	end

	% update
	Ax = Ax + step * Adir;
%	Cx = Cx + step * Cdir;
	x = x + step * ddir;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
	info(iter,:) = arg.userfun(x, iter, arg.userarg{:});

	% check norm(xnew-xold) / norm(xnew) vs threshold
	if arg.stop_diff_tol && ...
		norm_diff(step * ddir) / norm_diff(x) < arg.stop_diff_tol
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
function out = userfun_default(x, iter, varargin)
gamma = evalin('caller', 'gamma');
step = evalin('caller', 'step');
out = [gamma step cpu('etoc')];


function dot = dot_double(a, b)
dot = sum(a .* b, 'double'); % double accumulate


% pwls_pcg1_test
function pwls_pcg1_test
mask = true([8 7]); mask(1) = false;
A = 2 * Gdft('mask', mask); % orthogonal to permit analytical solution
xtrue = zeros(size(mask), 'single');
xtrue(end/2, round(end/2)) = 1;
y = A * xtrue(mask);
beta = 2^6;
% trick: using 'mat' not 'mex' option because 'mex' version has edge error
R = Reg1(mask, 'beta', beta, 'type_penal', 'mat');
hess = full(A' * A + R.C' * R.C);
xhat = hess \ (A' * y);
xhat = embed(xhat, mask);
im(xhat)

xinit = 0 * mask;
xpcg = pwls_pcg1(xinit(mask(:)), A, 1, y, R, 'niter', 100, ...
	'stop_grad_tol', 1e-6, 'stop_grad_norm', 2, 'chat', 1);
xpcg = embed(xpcg, mask);

im plc 2 2
im(1, xtrue)
im(2, xhat)
im(3, xpcg)
im(4, xpcg - xhat)
equivs(xpcg, xhat)
