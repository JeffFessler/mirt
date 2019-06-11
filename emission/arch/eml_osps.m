 function xs = eml_osps(x, Gb, yi, ci, ri, niter, pixmax, curv, relax0, chat)
%function xs = eml_osps(x, Gb, yi, ci, ri, niter, pixmax, curv, relax0, chat)
% E-ML-OSPS algorithm for emission Poisson problem
% (ordered subsets separable paraboloidal surrogates)
% model: Y_i ~ Poisson(c_i [G x]_i + r_i)
% in
%	x	[np,1]		initial estimate
%	Gb			Gblock object (see eml_osps_test.m)
%	yi,ci,ri [nb,na]	see em_fbp.m (for model too)
%	niter			# iterations
%	pixmax			upper constraint for pixel values
%	curv	'oc' for erdogan's optimal curvatures
%		'pc' for erdogan's fast precomputed curvatures,
%			which usually provides faster convergence,
%			but can be nonmonotone
%	relax0	[1] or [2]	relax0 or (relax0, relax_rate)
%	chat
% out
%	xs [np,niter]	updated image vectors each iteration
%
% Copyright Mar 2000, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

nblock = block_ob(Gb, 'n');
starts = subset_start(nblock);

if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end
if ~isvar('niter') || isempty(niter),	niter = 1;	end
if ~isvar('pixmax') || isempty(pixmax),	pixmax = inf;	end
if ~isvar('curv') || isempty(curv),	curv = 'oc';	end
if ~isvar('chat') || isempty(chat),	chat = true;	end

if ~isvar('relax0') || isempty(relax0)
	relax0 = 1;
end
if length(relax0) == 1
	relax_rate = 0;
elseif length(relax0) == 2
	relax_rate = relax0(2);
	relax0 = relax0(1);
else
	error relax
end

eml_check(yi, ci, ri);

[nb, na] = size(yi);

gi = sum(Gb')';		% g_i = sum_j g_ij
gi = reshape(gi, nb, na);

%
%	precomputed curvatures
%
denom = zeros(numel(x), nblock);
if streq(curv, 'pc')
	ni = eml_curvature(yi, ci, ri, [], [], curv);
	denom = Gb' * col(gi .* ni);
%	printf('ni range %g %g', min(ni(:)), max(ni(:)))
%	printf('denom range %g %g', min(denom(:)), max(denom(:)))
end


%
%	loop over iterations
%
xs = zeros(numel(x), niter);
x = max(x,0);
x = min(x,pixmax);
xs(:,1) = x;
for iter = 2:niter

	relax = 1;

	%
	% loop over subsets
	%
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Gb{iblock} * x;			% l=G*x "line integrals"
		li = reshape(li, nb, length(ia));
		yb = ci(:,ia) .* li + ri(:,ia);		% predicted meas. means

		% fix: need to be careful here with 0/0 -> 0
		dothi = ci(:,ia) .* (yi(:,ia) ./ yb - 1);

		% non-precomputed curvatures (notably, optimal curvature),
		% for ensured monotone increase
		if ~streq(curv, 'pc')
			ni = eml_curvature(yi(:,ia), ci(:,ia), ri(:,ia), li, yb, curv);
			denom = nblock * (Gb{iblock}' * col(gi(:,ia) .* ni));
		end

		grad = Gb{iblock}' * dothi(:);
		num = nblock * grad;
		x = x + relax * num ./ denom;	% relaxed update
		x = max(x,0);			% lower bound
		x = min(x,pixmax);		% upper bound
	end

	xs(:,iter) = x;
end
