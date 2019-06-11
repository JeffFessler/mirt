 function [xs, info] = qpwls_pcg1(x, G, W, yi, C, varargin)
%function [xs, info] = qpwls_pcg1(x, G, W, yi, C, [options])
%
% quadratic penalized weighted least squares via
% preconditioned conjugate gradients (PCG) algorithm
% cost(x) = (y-Gx)'W(y-Gx)/2 + x'C'Cx/2
% in
%	x	[np,1]		initial estimate
%	G	[nd,np]		system matrix
%	W	[nd,nd]		data weighting matrix
%	yi	[nd,1]		noisy data
%	C	[nc,np]		penalty 'derivatives' (R = \Half C'*C)
%
% options
%	niter			# total iterations (default: 1)
%	isave	[]		list of iterations to archive
%				(default: last iteration only)
%	userfun			user defined function handle (see default below)
%	precon	[np,np]		preconditioner (or object) (or 1)
%
% out
%	xs	[np,niter]	estimates each iteration
%	info	[niter, 3]	gamma, step, time
%
% Copyright Jan 1998, Jeff Fessler, The University of Michigan

if nargin < 5, help(mfilename), error(mfilename), end

% defaults
arg.precon = 1;
arg.niter = 1;
arg.isave = [];
arg.userfun = @userfun_default;
arg.key = 1;

arg = vararg_pair(arg, varargin);
if isempty(arg.isave), arg.isave = arg.niter; end
if streq(arg.isave, 'all'), arg.isave = 0:arg.niter; end

if ~isreal(yi) && ~isequal(arg.precon, 1)
	persistent warned
	if isempty(warned), warned = 0; end
	if ~warned
		warning 'not 100% sure about the complex preconditioned case'
		warned = 1;
	end
end

cpu etic

x = x(:);
np = length(x);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

%info = zeros(arg.niter, ?); % trick: do not initialize since size may change

%
% initialize projections
%
ticker(mfilename, 1, arg.niter)
Gx = G * x;
Cx = C * x;

%flops0 = flops;

%
% iterate
%
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	%
	% (negative) gradient
	%
	ngrad = G' * (W * (yi-Gx)) - C' * Cx;

	%
	% preconditioned gradient
	%
	pregrad = arg.precon * ngrad;

	%
	% direction
	%
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
			gamma = newinprod / oldinprod;	% Fletcher-Reeves
		end
		ddir = pregrad + gamma * ddir;
	end
	oldgrad = ngrad;
	oldinprod = newinprod;

	% check if descent direction
	if real(ddir' * ngrad) < 0
		warning 'wrong direction'
		if arg.key, keyboard, end
	end

	%
	% step size in search direction
	%
	Gdir = G * ddir;
	Cdir = C * ddir;

	denom = Gdir'*(W*Gdir) + Cdir'*Cdir;
	if denom == 0
		warning 'found exact solution??? step=0 now!?'
		step = 0;
	else
		denom = reale(denom, 'error', 'denom');
		step = (ddir' * ngrad) / denom;
%		step = reale(step, 'warn', 'step');
		step = real(step); % real step sizes seems only logical
	end
	if step < 0
		warning 'downhill?'
		if arg.key, keyboard, end
	end

	%
	% update
	%
	Gx	= Gx + step * Gdir;
	Cx	= Cx + step * Cdir;
	x	= x + step * ddir;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
        end
	info(1+iter,:) = arg.userfun(x);
end


% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default(x)
gamma = evalin('caller', 'gamma');
step = evalin('caller', 'step');
out = [gamma step cpu('etoc')];
