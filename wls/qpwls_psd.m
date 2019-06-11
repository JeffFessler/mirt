 function [xs, info] = qpwls_psd(x, A, W, yi, C, varargin)
%function [xs, info] = qpwls_psd(x, A, W, yi, C, [options])
%|
%| quadratic penalized weighted least squares (QPWLS) via
%| preconditioned steepest descent (PSD) algorithm
%| cost(x) = (y-Ax)'W(y-Ax) / 2 + x'C'Cx / 2
%|
%| in
%|	x	[np 1]		initial estimate
%|	A	[nd np]		system matrix
%|	W	[nd nd]		data weighting matrix
%|	yi	[nd 1]		noisy data
%|	C	[nc np]		penalty 'differencing matrix' (0 for unregularized)
%|
%| options
%|	niter			# total iterations (default: 1)
%|	isave	[]		list of iterations to archive
%|				(default: last iteration only)
%|	userfun	@		user defined function handle (see default below)
%|	userarg	{}		user arguments to userfun (default {})
%|	precon	[np np]		preconditioner (matrix or object) (or 1)
%|	dircheck 0|1		check descent direction? (default: 1)
%|				set to 0 to save time, if you dare...
%|	stop_threshold		stop iterations if norm(xnew-xold)/norm(xnew)
%|				is less than this unitless value.  default 0.
%|
%| out
%|	xs	[np niter]	estimates each iteration
%|	info	[niter 2]	step size, time each iteration
%|
%| Copyright Jun 2000, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), run_mfile_local('qpwls_psd_test_todo'), return, end
if nargin < 5, help(mfilename), error(mfilename), end

% defaults
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


cpu etic
if isempty(x), x = zeros(ncol(A),1); end

x = x(:);
np = length(x);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

%info = zeros(arg.niter, ?); % trick: do not initialize since size may change

% initialize projections
ticker(mfilename, 1, arg.niter)
Ax = A * x;
Cx = C * x;

% iterate
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	% (negative) gradient
	ngrad = A' * (W * (yi-Ax)) - C' * Cx;

	% preconditioned gradient: search direction
	ddir = arg.precon * ngrad;

	% check if descent direction
	if arg.dircheck && real(ddir' * ngrad) < 0
		warning 'wrong direction'
		if arg.key, keyboard, end
	end

	% step size in search direction
	Adir = A * ddir;
	Cdir = C * ddir;

	denom = Adir'*(W*Adir) + Cdir'*Cdir;
	if denom == 0
		warning 'found exact solution??? step=0 now!?'
		step = 0;
	else
		denom = reale(denom, 'error', 'denom');
		step = (ddir' * ngrad) / denom;
		step = real(step); % see a-complex.pdf
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
		xs(:, arg.isave == iter) = x;
	end
	info(1+iter,:) = arg.userfun(x, arg.userarg{:});

	% check norm(xnew-xold) / norm(xnew) vs threshold
	if arg.stop_threshold && ...
		norm(step * ddir) / norm(x) < arg.stop_threshold
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
step = evalin('caller', 'step');
out = [step cpu('etoc')];
