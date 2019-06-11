 function xs = wls_grpr(x, G, W, yi, D, xmin, xmax, niter)
%function xs = wls_grpr(x, G, W, yi, D, xmin, xmax, niter)
%
%	weighted least squares with constraint xmin <= x <= xmax
%	using gradient projection method (polyak:87 p. 207)
%		xnew = max(xmin, xold + D * \nabla J(xold))
%	cost function: J(x) = (y-Gx) W (y-Gx) / 2
%	in
%		x	[np,1]		initial estimate
%		G	[nd,np]		system matrix
%		W	[nd,nd]		weighting matrix
%		yi	[nd,1]		noisy measurements (e.g. sinogram)
%		D	[np,np]		preconditioning / step size matrix
%		xmin			minimum allowable x value
%		xmax			maximum allowable x value
%		niter			# of iterations
%	out
%		xs	[np,niter+1]	iterates
%
%	Copyright Dec. 2000,	Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

if ~isvar('W') || isempty(W)
	W = 1;
end
if ~isvar('niter') || isempty(niter)
	niter = 1;
end
if ~isvar('D') || isempty(D)
	D = 1;
end
if ~isvar('xmin') || isempty(xmin), xmin = 0; end
if ~isvar('xmax') || isempty(xmax), xmax = inf; end


% loop over iterations
xs = zeros(length(x), niter);
xs(:,1) = x;
for iter = 2:niter
%	printf('WLS-PG iteration %d', iter-1)

	if ~rem(iter,2)
		lgrad = G' * (W .* (yi - G * x));

		x = x + D * lgrad;	% the update!
	else

		x = max(x,xmin);	% lower bound
		x = min(x,xmax);	% upper bound
	end

	xs(:,iter) = x;
end
