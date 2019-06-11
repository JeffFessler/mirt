 function [xs, info] = pwls_pcg(x, G, W, yi, nder1, R, ...
		M, niter, stepper)
%function [xs, info] = pwls_pcg(x, G, W, yi, nder1, R, ...
%		M, niter, stepper)
%
% weighted least squares with convex non-quadratic penalty
% via preconditioned conjugate gradient algorithm
% cost(x) = (y-Gx)'W(y-Gx)/2 - n'(y-Gx) + R(x)
%
% in
%	x	[np,1]		initial estimate
%	G	[nd,np]		system matrix
%	W	[nd,nd]		data weighting matrix, usually diag_sp(wi)
%	yi	[nd,1]		noisy data
%	nder1	[nn,1]		linear term (obsolete: use "0")
%	R			penalty object (see Reg1.m)
%	M	[np,np]		preconditioner (use "1" if none)
%	niter			# total iterations
%	stepper			method for step-size line search
%				use {} for a good default
% out
%	xs	[np,niter]	estimates each iteration
%	info	[niter, 3]	gamma, step, time
%
% Copyright 1996-7, Jeff Fessler, The University of Michigan

warning 'pwls_pcg is obsolete; use pwls_pcg1 instead'

if nargin < 8, help(mfilename), error args, end

cpu etic

if isempty(nder1), nder1 = 0; end
if isempty(M), M = 1; end
if nder1 ~= 0, error 'nder1 = 0 required', end

if ~isvar('stepper') || isempty(stepper)
	stepper = {'qs', 3};	% quad surr with this # of subiterations
end

xs = zeros(length(x), niter);
xs(:,1) = x;

info = zeros(niter,3);

%
% initialize projections
%
Gx = G * x;

oldinprod = 0;

% iterate
for ii=2:niter
	ticker(mfilename, ii, niter)

	%
	% (negative) gradient
	%
	ngrad = G' * (W * (yi-Gx) - nder1);

	pgrad = R.cgrad(R, x);
	ngrad = ngrad - pgrad;

	%
	% preconditioned gradient
	pregrad = M * ngrad;

	% direction
	newinprod = ngrad' * pregrad;
	if ii == 2
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
	if real(ddir' * ngrad) < 0
		warning 'wrong direction'
		keyboard
%		ddir = pregrad;	% revert
%		oldinprod = 0;	% reset
	end

	% step size in search direction
	Gdir = G * ddir;
%	Cdir = R.C * ddir;

	% one step based on quadratic surrogate for penalty
	if streq(stepper{1}, 'qs1')
%		pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir); % cannot be?
		pdenom = (abs(ddir).^2)' * R.denom(R, x);
		denom = Gdir'*(W*Gdir) + pdenom;
		if denom == 0
			warning 'found exact solution???  step=0 now!?'
			step = 0;
		else
			step = real((ddir' * grad) / denom);
		end

	% iteratively minimize \Half || y-G (x+alf*ddir) ||_W^2 + R(x + alf*ddir)
	elseif streq(stepper{1}, 'qs')
		nsub = stepper{2};
		dGWGd = Gdir'*(W*Gdir);
		dGWr = Gdir'*(W*(yi-Gx));
		step = 0;
		for is=1:nsub
%			pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir);
			pdenom = (abs(ddir).^2)' * R.denom(R, x+step*ddir);
			denom = dGWGd + pdenom;
			pgrad = R.cgrad(R, x + step * ddir);
			step = step - (-dGWr + step * dGWGd + ddir' * pgrad) ...
				/ denom;
%			printf('%d-%d %g', ii, is, step)
		end

	else
		error 'bad stepper'
	end

	if step < 0
		warning('downhill?')
		keyboard
	end

	% update
	Gx	= Gx  + step * Gdir;
%	Cx	= Cx  + step * Cdir;
	x	= x + step * ddir;
	xs(:,ii) = x;

	info(ii,1) = gamma;
	info(ii,2) = step;
	info(ii,3) = cpu('etoc');	% accum. time
end
