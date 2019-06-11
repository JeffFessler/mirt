 function x = eql_emdp(x, G, yi, ci, ri, Asum, C)
%function x = eql_emdp(x, G, yi, ci, ri, Asum, C)
%	Alvaro De Pierro's modified EM algorithm
%	for quadratically penalized log-likelihood
%	image reconstruction from Poisson data
%	one iteration
%	Asum = A'1 = G' * c	(optional)
%	Penalty: R(x) = \sum_k ([Cx]_k)^2 / 2
%
%	Copyright 1998, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end
if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi(:)));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi(:)));
end
if ~isvar('Asum') || isempty(Asum)
	Asum = G' * ci;
end
if ~isvar('C') || isempty(C)
	C = 0;
end

eml_check(yi, ci, ri);
if any(x <= 0), error 'need x > 0', end

	yp = ci .* (G * x) + ri;		% predicted measurements
	if any(yi & ~yp), error 'model mismatch', end
	eterm = G' * (ci .* (yi ./ yp));

	%
	%	if \zkj = 1_{\ckj \neq 0} / \sum_j' 1_{\ckj' \neq 0}
	%	then \tilde{c}_j = \sum_k \ckj^2 / \zkj
	%
	%	in serial C we could be smarter and use a GEM approach
	% 
	n_per_k = sum(C' ~= 0)';
	if max(n_per_k) == 2
		% for consistency with ASPIRE using depierro=2 factor
%		warning 'using depierro = 2 version'
		Cfac = (C .* C)' * (2 * ones(size(n_per_k)));
	else
		Cfac = (C .* C)' * n_per_k;
	end

	x = eql_root(Cfac, (Asum + C'*(C*x) - Cfac .* x) / 2, x .* eterm);
