 function xs = tpl_os_sps(x, Gb, yi, bi, ri, R, niter, pixmax, curv, ...
		gi, denom, relax0, chat)
%function xs = tpl_os_sps(x, Gb, yi, bi, ri, R, niter, pixmax, curv, ...
%		gi, denom, relax0, chat)
%
% T-PL-OS-SPS algorithm for transmission Poisson problem
% (ordered subsets separable paraboloidal surrogates)
% in
%	x	[np,1]		initial estimate
%	Gb			Gblock object (see tpl_os_sps_test.m)
%	yi,bi,ri [nb,na]	see tr_fbp.m (for model too)
%	R		penalty structure (see Robject.m)
%				(or empty for ML case)
%	niter		# iterations
%	pixmax		upper constraint for pixel values or [lower, upper]
%			can be scalar (e.g. 'inf') or an array the size of x
%	curv		'oc' for erdogan's optimal curvatures
%			'pc' for erdogan's fast precomputed curvatures,
%				which usually provides faster convergence,
%				but can be nonmonotone
%	gi	[nd,1]	sum(G')'
%	denom	[np,1]		(if precomputed denominator)
%	relax0	[1] or [2]	relax0 or (relax0, relax_rate)
% out
%	xs [np,niter]	updated image vectors each iteration
%
% Copyright Mar 2000, Jeff Fessler, The University of Michigan
% 2002-2-20 update for Gblock object

% fix: add warning that this is supercede by pl_iot.m

if nargin < 3, help(mfilename), error(mfilename), end

Gb = block_ob(Gb, 'ensure'); % make it a block object (if not already)
nblock = block_ob(Gb, 'n');
starts = subset_start(nblock);

if ~isvar('bi') || isempty(bi)
	bi = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end

if ~isvar('R'), R = []; end
if ~isvar('niter')	|| isempty(niter),	niter = 2;	end
if ~isvar('pixmax')	|| isempty(pixmax),	pixmax = inf;	end
if length(pixmax) == 2
	pixmin = pixmax(1);
	pixmax = pixmax(2);
else
	pixmin = 0;
end
if ~isvar('curv')	|| isempty(curv),	curv = 'oc';	end
if ~isvar('chat')	|| isempty(chat),	chat = false;	end
if ~isvar('relax0')	|| isempty(relax0),	relax0 = 1;	end
if length(relax0) == 1
	relax_rate = 0;
elseif length(relax0) == 2
	relax_rate = relax0(2);
	relax0 = relax0(1);
else
	error relax
end

trl_check(yi, bi, ri);

if ~isvar('gi') || isempty(gi)
	gi = reshape(sum(Gb'), size(yi)); % g_i = sum_j g_ij
end

starts = subset_start(nblock);
[nb na] = size(yi);

% precompute denominator if needed
if (~isvar('denom') || isempty(denom)) && streq(curv, 'pc')
	ni = trl_curvature(yi, bi, ri, [], 'pc'); % precomputed curvatures

	%
	% one denominator shared by all subsets
	%
	denom = Gb' * col(gi .* ni);

	%
	% separate denominators for each subset
	% fix: this is inefficient, but consistent with 2d aspire!
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
% loop over iterations
%
xs = zeros(length(x), niter);
x = max(x,pixmin);
x = min(x,pixmax);
xs(:,1) = x;

for iter = 2:niter
	ticker(mfilename, iter, niter)

	relax = relax0 / (1 + relax_rate * (iter-2));

	%
	% loop over subsets
	%
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Gb{iblock} * x;	% l=G*x "line integrals"
		li = reshape(li, nb, length(ia));
		bel = bi(:,ia) .* exp(-li);
		yb = bel + ri(:,ia);	% predicted meas. means
		yb = yb + 100 * (yi(:,ia) == 0); % trick so that 0/0 -> 0.
		dothi = (1 - yi(:,ia) ./ yb) .* bel;

		% optimal curvatures (for ensured monotone increase)
		if ~streq(curv, 'pc')
			ni = trl_curvature(yi(:,ia), bi(:,ia), ri(:,ia), ...
				li, 'oc');
			denom = Gb{iblock}' * col(gi(:,ia) .* ni);
			denom = nblock * denom;
		end

		grad = Gb{iblock}' * dothi(:);

%%		for isub=1:1
			if isempty(R)
				num = nblock * grad;
				den = denom;
			else
				num = nblock * grad - R.cgrad(R, x);
				den = denom + R.denom(R, x);
			end
			x = x + relax * num ./ den;	% relaxed update
%%			x = x + ((x-xold) .* Lden + num) ./ den; % the update!
			x = max(x,pixmin);		% lower bound
			x = min(x,pixmax);		% upper bound
%%		end
	end

	if chat, printf('Range %g %g', min(x), max(x)), end
	xs(:,iter) = x;
end
