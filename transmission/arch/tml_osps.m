 function [xs, ni] = tml_osps(x, Gb, yi, bi, ri, niter, pixmax, curv, gi, denom)
%function [xs, ni] = tml_osps(x, Gb, yi, bi, ri, niter, pixmax, curv, gi, denom)
%	The T-ML-OSPS algorithm for transmission Poisson problem
%	(ordered subsets separable paraboloidal surrogates)
%	Input
%		Gb		Gblock object
%		see tr_fbp.m for model, G, yi, bi, ri
%		niter		# iterations
%		pixmax		upper constraint for pixel values
%			can be scalar (e.g. 'inf') or an array the size of x
%		curv	'oc' for erdogan's optimal curvatures
%			'pc' for erdogan's fast precomputed curvatures,
%				which usually provides faster convergence,
%				but can be nonmonotone
%		gi		sum(G')' (optional)
%		denom		[np,1] (if precomputed denominator)
%	Output
%		xs [np,niter]	updated image vectors each iteration
%
%	Copyright Mar 2000, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

if ~isvar('bi') || isempty(bi)
	bi = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end

if ~isa(Gb, 'Gblock'), error Gblock, end
nblock = block_ob(Gb, 'n');

trl_check(yi, bi, ri);

if ~isvar('niter') || isempty(niter)
	niter = 2;
end
if ~isvar('pixmax') || isempty(pixmax)
	pixmax = inf;
end
if ~isvar('curv') || isempty(curv)
	curv = 'oc';
end
if strcmp(curv, 'oc'),		is_oc = 1;
elseif strcmp(curv, 'pc'),	is_oc = 0;
else,	error 'curv unknown',	end

if ~isvar('gi') || isempty(gi)
	gi = reshape(sum(Gb'), size(yi)); % g_i = sum_j g_ij
end

starts = subset_start(nblock);
[nb, na] = size(yi);

% precompute denominator if needed
if (~isvar('denom') || isempty(denom)) && ~is_oc
	ni = trl_curvature_pre(yi, bi, ri); % precomputed curvatures

	% a single denominator shared by all subsets (scaled)
	if 1
		denom = Gb' * col(gi .* ni) / nblock;
	% separate denominators for each subset
	% fix: this is inefficient, but consistent with 2d aspire!
	else
		denom = zeros(numel(x), nblock);
		for iset=1:nblock
			iblock = starts(iset);
			ia = iblock:nblock:na;
			denom(:,iset) = Gb{iblock} * col(gi(:,ia) .* ni(:,ia));
		end
	end
end

pixmin = 0;

% loop over iterations
xs = zeros(numel(x), niter);
x = max(x,pixmin);
x = min(x,pixmax);
xs(:,1) = x;

for it=2:niter
%	printf('TML OSPS iteration %d', it)

	% loop over subsets
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		% l=G*x "line integrals"
		li = Gb{iblock} * x;

		li = reshape(li, nb, length(ia));
		yb = bi(:,ia) .* exp(-li) + ri(:,ia); % predicted meas. means
		dothi = bi(:,ia) .* (1 - yi(:,ia) ./ yb) .* exp(-li);

		if is_oc
			% optimal curvatures (for ensured monotone increase)
			ni = trl_curvature_opt(yi(:,ia), bi(:,ia), ri(:,ia), li, yb);
			denom = Gb{iblock}' * col(gi(:,ia) .* ni); % (:,iset)
		end

		grad = Gb{iblock}' * dothi(:);

%		x = x + grad ./ denom(:,iset);	% the update!
		x = x + grad ./ denom;		% the update!
		x = max(x,pixmin);		% enforce lower bound
		x = min(x,pixmax);		% enforce upper bound
	end

	xs(:,it) = x;
end
