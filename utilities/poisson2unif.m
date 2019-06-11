function [uu] = poisson2unif(count, lambda)

%	function [uu] = poisson2unif(count, lambda)
%	INPUT:
%		count: [Nx1]	measured counts
%		lambda: [Nx1]	hypothesized means
%	OUTPUT:
%		uu: [Nx1]	uniform numbers
%	see veklerov and llacer 1987

	[s0, d] = p_poisson(count, lambda);

	rand('uniform')
	uu = s0 + d .* rand(nrow(count), ncol(count));
