  function xs = eml_em(x, G, yi, ci, ri, Asum, niter)
%|function xs = eml_em(x, G, yi, ci, ri, Asum, niter)
%| E-ML-EM algorithm for image reconstruction from Poisson emission data
%| in
%|	x	[np 1]		initial image guess (column - see examples)
%|	Asum 			A'1 = G' * c	(optional) column sums
%|	see em_fbp.m for model, G, yi, ci, ri
%| out
%|	xs	[np niter]	image iterates
%|
%| Copyright 1998, Jeff Fessler, University of Michigan

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

if any(Asum <= 0), error 'Asum must be positive.  Adjust mask.', end

if ~isvar('niter') || isempty(niter)
	niter = 1;
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

	%
	% update:
	% x = x .* (A' * (y ./ (A * x + r))) ./ sum(A)'
	% where A = D(ci) G
	%

	yp = ci .* (G * x) + ri; % predicted measurements
	eterm = G' * (ci .* (yi ./ yp));
	x = x .* eterm ./ Asum;

	xs(:,iter) = x;
end
