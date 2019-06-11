 function [xs, ni] = tml_convex(x, G, yi, bi, ri, niter, pixmax, type)
%function [xs, ni] = tml_convex(x, G, yi, bi, ri, niter, pixmax, type)
% T-ML-Convex algorithms for transmission Poisson problem
% model: Y_i ~ Poisson(b_i exp(-[G x]_i) + r_i)
% in
%	x	[np,1]	initial guess
%	G	[nd,np]	transpose of system matrix
%	yi	[nd,1]	transmission sinogram
%	bi	[nd,1]	blank scan factors
%	ri	[nd,1]	background (randoms, scatter, crosstalk, etc)
%		bi,ri:	optional (can use empty matrices)
%		yi,bi,ri must have identical dimensions
%	niter		# of iterations
%	pixmax	upper constraint for pixel values
%		can be scalar (e.g. 'inf') or an array the size of x
%	type	'ps'	paraboloidal surrogate
%		'nr2'	newton raphson (2nd choice of zij's)
%		'log1'	-log method with amax
%		'log2'	-log with better choice fo zij's
% out
%	x [np,niter]	updated image vectors each iteration
%
% Copyright 2001-10-17, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

if ~isvar('bi') || isempty(bi)
	bi = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end
if ~isvar('niter') || isempty(niter)
	niter = 2;
end
if ~isvar('pixmax') || isempty(pixmax)
	pixmax = inf;
end
pixmin = 0;
if ~isvar('type') || isempty(type)
	type = 'ps';
end

trl_check(yi, bi, ri);

% column sums
gi = sum(G')';	% [nd,1] g_i = sum_j g_ij

if streq(type, 'log1')
	step = 1 / max(gi);
	printf('log1 step = %g', step)

elseif streq(type, 'log2')
	gj = sum(G)';				% [np,1] g_j = sum_i g_ij
	tmp = G * (1 ./ gj);			% [nd,1] sum_j g_ij / gj
	fj = max(tmp) * gj;
	tmp = G * (1 ./ fj);			% [nd,1] sum_j g_ij / gj
	tmp = full(tmp);
	printf('check inequality = %g, %g', min(tmp), max(tmp))
	step = 1 ./ fj;
	printf('log2 steps = %g, %g', min(step), max(step))

elseif streq(type, 'nr2')
	warning 'newton can be non-monotone'

elseif ~streq(type, 'ps')
	error 'unknown convex type'
end

% initialize
xs = zeros(length(x), niter);
x = max(x,pixmin);
x = min(x,pixmax);
xs(:,1) = x;

%
% loop over iterations
%
for it=2:niter
	li = reshape(G * x(:), size(yi));	% l=G*x "line integrals"
	bel = bi .* exp(-li);
	yb = bel + ri;			% predicted measurement means 

	if streq(type, 'log1') || streq(type, 'log2')
		top = G' * bel(:);
		Ni = bel ./ yb .* yi;
		bottom = G' * Ni(:);

		x = x + step .* log(top ./ bottom);	% the update!

	else
		dothi = (1 - yi ./ yb) .* bel;

		if streq(type, 'ps')
			Ni = yi .* bel ./ yb;			% \bar{N}_i
			ci = trl_curvature(Ni, bi, 0, li, 'oc');	% trick!
		elseif streq(type, 'nr2')
			ci = bel;
		else
			error bug
		end

                denom = G' * (gi .* ci(:));

		x = x + (G' * dothi(:)) ./ denom;	% the update!
	end

	x = max(x,pixmin);			% lower bound
	x = min(x,pixmax);			% upper bound

	xs(:,it) = x;
end
