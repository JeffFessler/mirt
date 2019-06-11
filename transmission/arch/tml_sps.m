 function [xs, ni] = tml_sps(x, Gt, yi, bi, ri, niter, pixmax, curv)
%function [xs, ni] = tml_sps(x, Gt, yi, bi, ri, niter, pixmax, curv)
%	One iteration of the ML-SPS algorithm for transmission Poisson problem
%	(separable paraboloidal surrogates)
%	model: Y_i ~ Poisson(b_i exp(-[G x]_i) + r_i)
%	Input
%		x	[np,1]	initial guess
%		Gt		transpose of system matrix
%		yi		transmission sinogram
%		bi		blank scan factors
%		ri		background (randoms, scatter, crosstalk, etc)
%		bi,ri:		optional (can use empty matrices)
%		yi,bi,ri	must have identical dimensions
%		pixmax		upper constraint for pixel values
%			can be scalar (e.g. 'inf') or an array the size of x
%		curv		'oc' for erdogan's optimal curvatures
%				'pc' for erdogan's fast precomputed curvatures,
%					which usually gives faster convergence,
%					but can be nonmonotone
%				'nr' newton curvatures, can be nonmonotone
%	Output
%		x [np,niter]	updated image vectors each iteration
%
%	fix: the really slick way to do this would be to use
%	the fast precomputed denominator and just backtrack to the monotone
%	version on those rare occasions when it goes downhill
%
%	Copyright 2000-3-01	Jeff Fessler	The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

[nb, na] = size(yi);

if (nargin < 4 || isempty(bi))
	bi = ones(size(yi));
end
if (nargin < 5 || isempty(ri))
	ri = zeros(size(yi));
end
if (nargin < 6 || isempty(niter))
	niter = 2;
end
if (nargin < 7 || isempty(pixmax))
	pixmax = inf;
end
if (nargin < 8 || isempty(curv))
	curv = 'oc';
end

trl_check(yi, bi, ri);

	gi = sum(Gt)';	% g_i = sum_j g_ij

	if strcmp(curv, 'pc')
		ni = trl_curvature_pre(yi, bi, ri);	% precomputed
		denom = Gt * (gi .* ni(:));
	elseif strcmp(curv, 'nr')
		warning 'newton curvatures can be non-monotone'
	elseif ~strcmp(curv, 'oc')
		error 'curv not implemented'
	end

xs = zeros(length(x), niter);
x = max(x,0);
x = min(x,pixmax);
xs(:,1) = x;

%
%	loop over iterations
%
for ii=2:niter
	li = reshape(Gt' * x, size(yi));	% l=G*x "line integrals"
	bel = bi .* exp(-li);
	yb = bel + ri;			% predicted measurement means 

	dothi = (1 - yi ./ yb) .* bel;

	if strcmp(curv, 'oc')
		%	optimal curvatures (for ensured monotone increase)
		ni = trl_curvature(yi, bi, ri, li, 'oc');
		denom = Gt * (gi .* ni(:));
	elseif strcmp(curv, 'nr')
		ni = (1 - ri.*yi./yb.^2) .* bel;
		denom = Gt * (gi .* ni(:));
	end

	x = x + (Gt * dothi(:)) ./ denom;	% there's the update!
	x = max(x,0);				% enforce nonnegativity
	x = min(x,pixmax);			% enforce upper bound constraint

	xs(:,ii) = x;
end
