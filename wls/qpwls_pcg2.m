 function [xs, info] = qpwls_pcg2(x, T, bb, C, varargin)
%function [xs, info] = qpwls_pcg2(x, T, bb, C, [options])
%|
%| quadratic penalized weighted least squares via
%| preconditioned conjugate gradients (PCG) algorithm
%| cost: Psi(x) = (y-Ax)'W(y-Ax)/2 + x'C'Cx/2
%|		= -x'b + x'Tx/2 + x'C'Cx/2, where T = A'WA, b = A'W*y
%|
%| Essentially this routine solves (T + C'C) x = b
%| assuming that the matrix T is hermitian symmetric.
%|
%| in
%|	x	[np 1]		estimate
%|	T	[nd np]		system matrix
%|	bb	[np 1]		A'*W*y ("backprojection vector")
%|	C	[nc np]		penalty 'derivatives' (R = \Half C'*C)
%|
%| options
%|	niter			# total iterations (default: 1)
%|	isave	[]		list of iterations to archive
%|				(default: last iteration only)
%|	userfun	@		user defined function handle (see default below)
%}	userarg	 {}		user arguments to userfun (default {})
%|					taking arguments (x, arg.userarg{:})
%|	precon	[np np]		preconditioner (or object) (or 1)
%|	dircheck 0|1		check descent direction? (default: 1)
%|				set to 0 to save time, if you dare...
%|	stop_threshold		stop iterations if norm(xnew-xold)/norm(xnew)
%|				is less than this unitless value.  default 0.
%|
%| out
%|	xs	[np niter]	estimates each iteration
%|	info	[niter ?]	gamma, step, time
%|
%| Copyright Jan 1998, Jeff Fessler & Hugo Shi, University of Michigan

if nargin == 1 && streq(x, 'test')
	run_mfile_local('qpwls_pcg2_test'), return
end
if nargin < 3, help(mfilename), error(mfilename), end

if ~isvar('C') || isempty(C), C = 0; end

arg.precon = 1;
arg.niter = 1;
arg.isave = [];
arg.userfun = @userfun_default;
arg.userarg = {};
arg.key = 1;
arg.dircheck = true; % default is to check descent direction
arg.stop_threshold = 0;

arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);

if ~isreal(bb) && ~isequal(arg.precon, 1)
	persistent warned
	if isempty(warned), warned = 0; end
	if ~warned
		warning 'i am not sure if the complex case is 100% correct'
		warned = 1;
	end
end

cpu etic

%
% initialize
%
x = x(:);
np = length(x);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

%info = zeros(niter,?);

Cx = C * x;
Tx = T * x;

%flops0 = flops;


% iterate
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	% (negative) gradient
	ngrad = bb - Tx - C' * Cx;

	% preconditioned gradient
	pregrad = arg.precon * ngrad;

	% direction
	newinprod = ngrad' * pregrad;
	% fix: should i take the real part?
	newinprod = reale(newinprod, 'warn', 'inprod');
	if iter == 1
		ddir = pregrad;
		gamma = 0;
	else
	%	gamma = (newinprod - oldgrad' * pregrad) / oldinprod;
		if (oldinprod == 0)
			warning 'inprod=0. going nowhere!'
			gamma = 0;
		else
			gamma = newinprod / oldinprod; % Fletcher-Reeves
		end
		ddir = pregrad + gamma * ddir;
	end
	oldgrad = ngrad;
	oldinprod = newinprod;

	% check if descent direction
	if arg.dircheck && real(ddir' * ngrad) < 0
		warning 'wrong direction'
		if arg.key, keyboard, end
	end

	% step size in search direction
	Cdir = C * ddir;
	Tdir = T * ddir;

	denom = ddir' * Tdir + Cdir' * Cdir;
	if denom == 0
		warning 'found exact solution??? step=0 now!?'
		step = 0;
	else
		step = (ddir' * ngrad) / denom;
		% fix: should i take the real part?
%		step = reale(step, 'warn', 'step');
		step = real(step); % real step seems only natural
	end
	if step < 0
		warning 'downhill?' % todo: fix to not warn if not descent
		if arg.key, keyboard, end
	end

	% update
	Tx = Tx + step * Tdir;
	Cx = Cx + step * Cdir;
	x = x + step * ddir;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
        info(1+iter,:) = arg.userfun(x, arg.userarg{:});

	% check norm(xnew-xold) / norm(xnew) vs threshold
	if arg.stop_threshold && ...
		norm(step * ddir) < norm(x) * arg.stop_threshold
		if isequal(arg.isave, arg.niter) % saving last iterate only?
			xs(:,1) = x; % be sure to save this 'final' iterate
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
