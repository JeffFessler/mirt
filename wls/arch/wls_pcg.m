 function [xs, info] = wls_pcg(G, W, yy, x0, M, niter, xmin, mask)
%function [xs, info] = wls_pcg(G, W, yy, x0, M, niter, xmin, mask)
%
% quadratic penalized weighted least squares via
% preconditioned conjugate gradients algorithm
% cost:	Psi(x) = (y-Gx)'W(y-Gx)/2
% in
%	G	[nd,np]		system matrix
%	W	[nd,nd]		data weighting matrix
%	yy	[nd,1]		noisy data
%	x0	[np,1] or [nx,ny]	initial estimate
%	M	[np,np]		preconditioner (or object) (or 1)
%	niter			# total iterations
%	xmin			minimum pixel value
%				DOES NOT WORK! JUST FOR DEMONSTRATION!
%	mask	[nx,ny]		which pixels are updated
% out
%	xs	[np,niter]	estimates each iteration
%	info	[niter, 3]	gamma, step, time
%
% Copyright 1998-07-03, Jeff Fessler, The University of Michigan

if nargin < 6 || nargin > 8, help(mfilename), error(mfilename), end

if isvar('mask') && ~isempty(mask), chat = 1; else, chat = 0; end
if ~isvar('xmin') || isempty(xmin)
	xmin = -Inf;
else
	warning('using xmin ruins convergence!')
end

np = ncol(G);
if numel(x0) ~= np && isvar('mask')
	x0 = x0(mask(:));
end

yy = yy(:);
xs = zeros(np, niter);
xs(:,1) = x0;

if chat
	clf, im(121, embed(x0,mask), 'x0')
end

info	= zeros(niter,3);

%
% initialize projections
%
Gx = G * xs(:,1);

tic

%
% iterate
%
oldinprod = 0;
ii = 2;
while ii <= niter
	old = xs(:,ii-1);

	%
	% (negative) gradient
	%
	ngrad = G' * (W * (yy-Gx));

	%
	% preconditioned gradient
	%
	pregrad = M * ngrad;

	%
	% direction
	%
	newinprod = ngrad' * pregrad;
	if oldinprod == 0
		ddir = pregrad;
		gamma = 0;
	else
%		gamma = (newinprod - oldgrad' * pregrad) / oldinprod;
		gamma = newinprod / oldinprod;
		ddir = pregrad + gamma * ddir;
	end
%	oldgrad = ngrad;
	oldinprod = newinprod;

	Gdir = G * ddir;

	% check if descent direction
	if ddir' * ngrad < 0
		warning('wrong direction')
		keyboard
	end

	%
	% step size in search direction
	%
	step_denom = Gdir'*(W*Gdir);
	if step_denom == 0
		warning 'zero denominator for step?'
		keyboard
	end
	step = (ddir' * ngrad) / step_denom;
	if step < 0
		warning('downhill?')
		keyboard
	end

	%
	% update
	%
	Gx	= Gx + step * Gdir;
	new	= old + step * ddir;
	xs(:,ii) = new;

	if chat && (ii < 10 || rem(ii,10)==0)
		im(121, embed(new,mask), 'xnew')
		title(sprintf('%d step=%g gamma=%g', ii, step, gamma))
		if isobject(M)
			ptype = M.type;
		else
			ptype = 'unknown';
		end
		im(122, embed(pregrad,mask), ['PreGrad ' ptype])
		drawnow
	end

	info(ii,1) = gamma;
	info(ii,2) = step;
	info(ii,3) = toc;		% accum. time

	if isinf(xmin)
		ii = ii + 1;
	else
		xs(:,ii+1) = max(xs(:,ii), xmin);	% enforce constraint
		Gx = G * xs(:,ii+1);			% serious extra work!!
		ii = ii + 2;
	end
end
