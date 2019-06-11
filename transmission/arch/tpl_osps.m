 function xs = tpl_osps(x, Gb, yi, bi, ri, P, niter, pixmax, curv, gi, denom)
%function xs = tpl_osps(x, Gb, yi, bi, ri, P, niter, pixmax, curv, gi, denom)
%	The T-PL-OSPS algorithm for transmission Poisson problem
%	(ordered subsets separable paraboloidal surrogates)
%	Input
%		x	[np,1]	initial guess
%		Gb		Gblock system object
%		see tr_fbp.m for model, G, yi, bi, ri
%		P		penalty structure: P.C,wpot
%		niter		# iterations
%		pixmax		upper constraint for pixel values
%			can be scalar (e.g. 'inf') or an array the size of x
%		curv	'oc' for erdogan's optimal curvatures
%			'pc' for erdogan's fast precomputed curvatures,
%				which usually provides faster convergence,
%				but can be nonmonotone
%		gi		sum(G')'
%		denom		[np,1] (if precomputed denominator)
%	Output
%		xs [np,niter]	updated image vectors each iteration
%
%	Copyright Mar 2000, Jeff Fessler, The University of Michigan
%	2002-2-20 update for Gblock object

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

if ~isvar('P') || isempty(P)
	error 'P required'
end
if ~isvar('niter') || isempty(niter)
	niter = 2;
end
if ~isvar('pixmax') || isempty(pixmax)
	pixmax = inf;
end
pixmin = 0;
if ~isvar('curv') || isempty(curv)
	curv = 'oc';
end
if strcmp(curv, 'oc'),		is_oc = 1;
elseif strcmp(curv, 'pc'),	is_oc = 0;
else,	error 'curv unknown',	end

if ~isvar('gi') || isempty(gi)
	gi = reshape(sum(Gb'), size(yi));	% g_i = sum_j g_ij
end

starts = subset_start(nblock);
[nb na] = size(yi);

%	precompute denominator if needed
if (~isvar('denom') || isempty(denom)) && ~is_oc
	ni = trl_curvature(yi, bi, ri, [], 'pc');	% precomputed curvatures

	%
	%	a single denominator shared by all subsets (scaled)
	%
	denom = Gb' * col(gi .* ni) / nblock;

	%
	%	separate denominators for each subset
	%	fix: this is inefficient, but consistent with 2d aspire!
	%
	if 0
		denom = zeros(numel(x), nblock);
		for iset=1:nblock
			iblock = starts(iset);
			ia = iblock:nblock:na;
			denom(:,iset) = Gb{iblock}' * col(gi(:,ia) .* ni(:,ia));
		end
	end
end

%
%	if \zkj = 1_{\ckj \neq 0} / \sum_j' 1_{\ckj' \neq 0}
%	then \tilde{c}_j = \sum_k \ckj^2 / \zkj * wpot(0)	<- max curvature
%
%	in serial C we could be smarter and use a GEM approach
%
C = P.C;
%n_per_k = sum(C' ~= 0)';
%if max(n_per_k) == 2
%	% for consistency with ASPIRE using depierro=2 factor
%%	warning 'using depierro = 2 version'
%	n_per_k(:) = 2;
%end
%	Cfac = (C .* C)' * (n_per_k .* P.wpot(0));	% max curvature
CCnt = (spdiag(sum(P.C' ~= 0)) * (P.C .^2))';

%
%	loop over iterations
%
xs = zeros(length(x), niter);
x = max(x,pixmin);
x = min(x,pixmax);
xs(:,1) = x;

for it=2:niter
%	printf('TPL OSPS iteration %d', it)

	%
	%	loop over subsets
	%
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Gb{iblock} * x;	% l=G*x "line integrals"
		li = reshape(li, nb, length(ia));
		bel = bi(:,ia) .* exp(-li);
		yb = bel + ri(:,ia);	% predicted meas. means
		dothi = (1 - yi(:,ia) ./ yb) .* bel;

		if is_oc
			% optimal curvatures (for ensured monotone increase)
			ni = trl_curvature(yi(:,ia), bi(:,ia), ri(:,ia), ...
				li, 'oc');
			denom = Gb{iblock}' * col(gi(:,ia) .* ni);
		end

		grad = Gb{iblock}' * dothi(:);

%%		for isub=1:1
			Cx = C * x;
			wx = P.wpot(P.wt, Cx);
			pgrad = C' * (wx .* Cx);
			num = nblock * grad - pgrad;

%			Cfac = (C .* C)' * (n_per_k .* P.wpot(Cx));
			Cfac = CCnt * wx;

%			den = nblock * denom(:,iset) + Cfac;
			den = nblock * denom + Cfac;
			x = x + num ./ den;		% the update!
%%			x = x + ((x-xold) .* Lden + num) ./ den; % the update!
			x = max(x,pixmin);		% lower bound
			x = min(x,pixmax);		% upper bound
%%		end
	end

	xs(:,it) = x;
end
