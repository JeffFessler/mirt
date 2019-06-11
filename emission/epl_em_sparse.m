 function xs = epl_em_sparse(x, G, yi, ci, ri, Asum, niter, Beta, gammaj, norm)
%function xs = epl_em_sparse(x, G, yi, ci, ri, Asum, niter)
% E-PL-EM algorithm for image reconstruction from Poisson emission data,
% regularized by the ell-<norm> norm, encouraging sparse solutions.
%
% in
%	x	[np,1]		initial image guess (column - see examples)
%	Asum			A'1 = G' * c	(optional) column sums
%	see em_fbp.m for model, G, yi, ci, ri
%	gammaj - vector of gammaj parameters
%	norm - 0 or 1 for ell-0 or ell-1 reconstruction
%
% out
%	xs	[np,niter]	image iterates
%
% Copyright 1998-2008, Jeff Fessler & Dan Lingenfelter, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

if ~isvar('ci') || isempty(ci)
	ci = 1;	% ones(size(yi(:)));
end
if ~isvar('ri') || isempty(ri)
	ri = 0; % zeros(size(yi(:)));
end
if ~isvar('Asum') || isempty(Asum)
	if length(ci) == 1
		Asum = sum(G)' * ci;
	else
		Asum = G' * ci;
	end
end

if ~isvar('norm') || isempty(norm)
	norm = 1;
end

if ~isvar('gammaj') || isempty(gammaj)
	gammaj = zeros(size(x));
end

if any(Asum <= 0), error 'Asum must be positive.  Adjust mask.', end
if and(norm ~= 0, norm ~= 1), error 'Norm must be either 0 or 1', end

if ~isvar('niter') || isempty(niter)
	niter = 1;
end

if ~isvar('Beta') || isempty(Beta)
	Beta = 0;
end

eml_check(yi, ci, ri);

if any(x <= 0), error 'need x > 0', end

xs = zeros(numel(x), niter);
xs(:,1) = x(:);

%
% loop over iterations
%
for iter = 2:niter
	ticker(mfilename, iter, niter)

	yp = ci .* (G * x) + ri;	% predicted measurements
	eterm = G' * (ci .* (yi ./ yp));

	if norm == 1
		x = (x + gammaj) .* eterm ./ (Asum + Beta) - gammaj;
	else
		xhat = (x + gammaj) .* eterm ./ (Asum) - gammaj;
		x = xhat .* (eterm .* (x + gammaj) .* ...
			log((xhat + gammaj) ./ gammaj) - Asum.*xhat > Beta);
	end

	xs(:,iter) = x;
end
