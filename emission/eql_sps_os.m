 function xs = eql_sps_os(x, Gb, yi, ci, ri, R, ...
		niter, pixmax, curv, relax0, chat)
%function xs = eql_sps_os(x, Gb, yi, ci, ri, R, ...
%		niter, pixmax, curv, relax0, chat)
%
% E-QPL-SPS-OS algorithm for emission Poisson problem
% (ordered subsets separable paraboloidal surrogates)
% in
%	x	[np,1]		initial estimate
%	Gb	[nd,np]		Gblock object (see eql_osps_test.m)
%	yi,ci,ri [nb,na]	see em_fbp.m (for model too)
%	R			penalty structure (see Robject.m)
%					(or empty for ML case)
%	niter			# iterations
%	pixmax			upper constraint for pixel values
%	curv		'oc' for erdogan's optimal curvatures
%			'pc' for erdogan's fast precomputed curvatures,
%				which usually provides faster convergence,
%				but can be nonmonotone.
%	relax0	[1] or [2]	relax0 or (relax0, relax_rate)
% out
%	xs [np,niter]	updated image vectors each iteration
%
% Copyright Mar 2000, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

Gb = block_ob(Gb, 'ensure'); % make it a block object (if not already)
nblock = block_ob(Gb, 'n');
starts = subset_start(nblock);

if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end

if ~isvar('R'), R = []; end
if ~isvar('niter')	| isempty(niter),	niter = 1;	end
if ~isvar('pixmax')	| isempty(pixmax),	pixmax = inf;	end
if ~isvar('curv')	| isempty(curv),	curv = 'oc';	end
if ~isvar('chat')	| isempty(chat),	chat = false;	end

if ~isvar('relax0')	| isempty(relax0),	relax0 = 1;	end
if length(relax0) == 1
	relax_rate = 0;
elseif length(relax0) == 2
	relax_rate = relax0(2);
	relax0 = relax0(1);
else
	error relax
end

eml_check(yi, ci, ri, 'os', nblock);

[nb, na] = size(yi);

gi = sum(Gb')';		% g_i = sum_j g_ij
gi = reshape(gi, nb, na);

%
% precomputed curvatures
%
if streq(curv, 'pc')
	ni = ci.^2 ./ max(yi,1);	% precomputed
	ni = eml_curvature(yi, ci, ri, [], [], 'pc');

	% efficient single denominator consistent with aspire
	denom = Gb' * col(gi .* ni);
elseif ~streq(curv, 'oc')
	error 'curv not implemented'
end


%
% loop over iterations
%
xs = zeros(numel(x), niter);
x = max(x,0);
x = min(x,pixmax);
xs(:,1) = x;

for iter = 2:niter
	if chat, printf('E-QL-SPS-OS iteration %d', iter-1), end

	relax = relax0 / (1 + relax_rate * (iter-2));

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
		grad = Gb{iblock}' * dothi(:);

		if streq(curv, 'oc')
			% optimal curvatures (for monotone increase)
			ni = eml_curvature(yi(:,ia), ci(:,ia), ri(:,ia), ...
					li, yb, 'oc');
			denom = Gb{iblock}' * col(gi(:,ia) .* ni);
			denom = nblock * denom;
		end

		if isempty(R)
			num = nblock * grad;
			den = denom;
		else
			num = nblock * grad - R.cgrad(R, x);
			den = denom + R.denom(R, x);
		end

		x = x + relax * num ./ den;	% relaxed update
		x = max(x,0);			% lower bound
		x = min(x,pixmax);		% upper bound
	end

	if chat, printf('Range %g %g', min(x), max(x)), end
	xs(:,iter) = x;
end
