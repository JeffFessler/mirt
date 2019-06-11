 function [xs, params, times] = pwls_qs_pcg(G, W, yy, x0, R, M, niter)
% ?
% weighted least squares with convex non-quadratic penalty
% preconditioned conjugate gradients algorithm for quadratic surrogates
%
% THIS IS NEW AND UNTESTED!!!
%
% cost(x) = (y-Gx)'W(y-Gx)/2 + R(x)
% in
%	G	[nd,np]		system matrix
%	W	[nd,nd]		data weighting matrix, usually diag_sp(wi)
%	yy	[nd,1]		noisy data
%	x0	[np,1] or [nx,ny]	initial estimate
%	R			roughness penalty object, see Reg1.m
%	M	[np,np]		preconditioner
%	niter			# total iterations
%?	stepper			method for step-size line search
%?	nsub			# subiterations for line search
% out
%	xs	[np,niter+1]	estimates each iteration
%?	params	[niter, 2]	gamma, step
%?	times	[niter, 2]	iteration time, flops
%
% Copyright 2003-4-24, Jeff Fessler, The University of Michigan

if nargin < 4, help mfilename, error mfilename, end

if ~isvar('stepper'), stepper = 'hq0'; end
if ~isvar('nsub'), nsub = 1; end

if ~isvar('R') || isempty(R)
	Rgrad = 0;	% unregularized default
	Rwpot = 0;
end

np = ncol(G);
yy = yy(:);
xs = zeros(length(x0), niter+1);
xs(:,1) = x0;

%	params	= zeros(niter,2);
%	times	= zeros(niter,2);

%
% initialize projections
%
Gx = G * xs(:,1);
Cx = C * xs(:,1);

tic
% flops0 = flops;

%
% iterate
%
for ii=1:niter
	old = xs(:,ii);

	%
	% (negative) gradient
	%
	grad = G' * (W * (yy-Gx));

	if ~isempty(R)
		Rgrad = R.cgrad(R, x);
		Rwpot = R.wpot(R, Cx);
	end

	grad = grad - Rgrad;

	%
	% preconditioned gradient
	%
	if isempty(M)
		pregrad = grad;
	else
		pregrad = M * grad;
	end

	%
	% direction
	%
	newinprod = grad' * pregrad;
	if ii == 1
		ddir = pregrad;
		gamma = 0;
	else
		if 1
			gamma = 0;	% HARDWIRE TO PSD FOR NOW!
%		elseif isquad
%			gamma = newinprod / oldinprod;	% quadratic
		else
			gamma = (newinprod - oldgrad' * pregrad) / oldinprod;
		end
		ddir = pregrad + gamma * ddir;
	end
	oldgrad = grad;
	oldinprod = newinprod;

	Gdir = G * ddir;
	Cdir = C * ddir;

	% check if descent direction
	if ddir' * grad < 0
		warning('wrong direction')
		keyboard
	end

	%
	% step size in search direction, by minimizing quadratic surrogate
	%
	step = (ddir' * grad) / (Gdir'*(W*Gdir) + Cdir'*(Rwpot .* Cdir));
	if step < 0
		warning('downhill?')
		keyboard
	end

	%
	% update
	%
	Gx	= Gx  + step * Gdir;
	Cx	= Cx  + step * Cdir;
	new	= old + step * ddir;
	xs(:,ii+1) = new;

%	params(ii,1) = gamma;
%	params(ii,2) = step;
%	times(ii,1) = toc;		% accum. time
%	times(ii,2) = flops;		% accum. flops

end

% times(:,2) = times(:,2) - flops0;
