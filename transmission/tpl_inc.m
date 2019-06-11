 function xs = tpl_inc(x, Gb, yi, bi, ri, R, niter, pixmax, curv, ...
		subout, enhance, gi, denom, relax0, chat)
%function xs = tpl_inc(x, Gb, yi, bi, ri, R, niter, pixmax, curv, ...
%		subout, enhance, gi, denom, relax0, chat)
%
% TRIOT: incremental SPS algorithm for transmission Poisson problem
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
%			'mc' for maximum curvature
%			'pc' for erdogan's fast precomputed curvatures,
%				which usually provides faster convergence,
%				but can be nonmonotone.
%	subout	[1]	nonzero value: return all subiterates
%			otherwise, return only iterates
%	enhance [1]	nonzero value: Hsiao's enhancement
%			otherwise, no enhancement
%	gi	[nd,1]	sum(G')'
%	denom	[np,1]		(if precomputed denominator)
%	relax0	[1] or [2]	relax0 or (relax0, relax_rate)
% out
%	xs [np,niter]	updated image vectors each iteration
%
% Copyright Mar 2000, Jeff Fessler, The University of Michigan
% 2002-2-20 update for Gblock object
% Modified for incremental version by Sangtae Ahn, from ~sangtaea/trans.

epsilon = 1e-10;

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
if ~isvar('subout')	|| isempty(subout),	subout = 0;	end
if ~isvar('enhance')	|| isempty(enhance),	enhance = 0;	end
if ~isvar('chat')	|| isempty(chat),	chat = false;	end
if ~isvar('relax0')	|| isempty(relax0),	relax0 = 1;	end

if streq(curv, 'mc')
	subscale = @(iset) iset;
else
	subscale = @(iset) 1;
end

if length(relax0) == 1
	relax_rate = 0;
elseif length(relax0) == 2
	relax_rate = relax0(2);
	relax0 = relax0(1);
else
	error relax
end

if length(niter) == 2
	osspsiter = niter(1);
	niter = niter(2);
else
	osspsiter = 1;
end

trl_check(yi, bi, ri);

if ~isvar('gi') || isempty(gi)
	gi = reshape(sum(Gb'), size(yi));	% g_i = sum_j g_ij
end

starts = subset_start(nblock);
[nb na] = size(yi);

% precompute denominator if needed
if (~isvar('denom') || isempty(denom)) && (streq(curv, 'pc') || streq(curv, 'mc'))
	ni = my_trl_curvature(yi, bi, ri, [], curv); % precomputed curvatures
	%
	% separate denominators for each subset
	%
	if streq(curv, 'mc')
		denom = zeros(numel(x), nblock);
		for iset=1:nblock
			iblock = starts(iset);
			ia = iblock:nblock:na;
			denom(:,iset) = Gb{iblock}' * col(gi(:,ia) .* ni(:,ia));
		end
	else
	%
	% one denominator shared by all subsets
	%
		denom = Gb' * col(gi .* ni) / nblock;
	end
end

% "precomputed curvature"
if ~streq(curv, 'pc')
	ni = my_trl_curvature(yi, bi, ri, [], 'pc');
	denom_pc = Gb' * col(gi .* ni) / nblock;
else
	denom_pc = denom;
end


if subout
	xs = zeros(length(x), (niter-1)*nblock+1);
else
	xs = zeros(length(x), niter);
end

x = max(x,pixmin);
x = min(x,pixmax);
xs(:,1) = x;
if subout, idx = 2; end
num_store = zeros(length(x), nblock);
den_store = zeros(length(x), nblock);

%
% initialization by os-sps-pc
%
for iter = 1:osspsiter
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Gb{iblock} * x;	% l=G*x "line integrals"
		li = reshape(li, nb, length(ia));
		bel = bi(:,ia) .* exp(-li);
		yb = bel + ri(:,ia);	% predicted meas. means
		dothi = (1 - yi(:,ia) ./ yb) .* bel;

		if streq(curv, 'oc')
			ni = my_trl_curvature(yi(:,ia), bi(:,ia), ri(:,ia), ...
				li, 'oc');
			denom = Gb{iblock}' * col(gi(:,ia) .* ni);
		end

		if isempty(R)
			grad = Gb{iblock}' * dothi(:);
			den_osps = max(denom_pc, epsilon);
			den_store(:, iset) = max(denom(:,subscale(iset)), ...
						epsilon);
		else
			grad = Gb{iblock}' * dothi(:) - R.cgrad(R, x)/nblock;
			den_osps = max(denom_pc + R.denom(R, x)/nblock, ...
					epsilon);
			den_store(:, iset) = max(denom(:,subscale(iset)) + ...
					R.denom(R, x)/nblock, epsilon);
		end

		if iter == osspsiter
			num_store(:, iset) = x.*den_store(:, iset) + grad;
		end
		x = x + grad ./ den_osps;
		x = max(x,pixmin);
		x = min(x,pixmax);
		if subout
			xs(:,idx) = x;
			idx = idx + 1;
		end
	end
	if ~subout, xs(:,iter+1) = x;	end
end
num = sum(num_store, 2);
den = sum(den_store, 2);

if 1,		% optional
	x = num ./ den;
	x = max(x,pixmin); x = min(x,pixmax);
	if subout
		xs(:,idx-1) = x;
	else
		xs(:,iter+1) = x;
	end
end


alphidx = 1;
for iter = (osspsiter+2):niter
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
		dothi = (1 - yi(:,ia) ./ yb) .* bel;

		% optimal curvatures (for ensured monotone increase)
		if streq(curv, 'oc')
			ni = my_trl_curvature(yi(:,ia), bi(:,ia), ri(:,ia), ...
				li, 'oc');
			denom = Gb{iblock}' * col(gi(:,ia) .* ni);
		end

		num = num - num_store(:, iset);
		den = den - den_store(:, iset);

		if isempty(R)
			grad = Gb{iblock}' * dothi(:);
			den_store(:, iset) = max(denom(:,subscale(iset)), ...
						epsilon);
		else
			grad = Gb{iblock}' * dothi(:) - R.cgrad(R, x)/nblock;
			den_store(:, iset) = max(denom(:,subscale(iset)) + ...
						R.denom(R, x)/nblock, epsilon);
		end

		num_store(:, iset) = x.*den_store(:, iset) + grad;

		num = num + num_store(:, iset);
		den = den + den_store(:, iset);

		xold = x;

		x = num ./ den;
		x = max(x,pixmin); x = min(x,pixmax);

		if enhance
			xos = xold + grad ./ (den_store(:,iset) - ...
				denom(:,subscale(iset)) + denom_pc);
			xos = max(xos,pixmin);
			xos = min(xos,pixmax);

			c2 = - sum((xos - x).^2 .* den) / 2;
			c1 = num'*(xos - x) - sum(x.*(xos-x).*den);
			c0 = num'*(x - xold) - sum(x.^2.*den - xold.^2.*den)/2;

			if 1
				if c0+c1+c2 > 0
					alpha = 1;
				else
					alpha = eql_root(-c2, -c1/2, 0.95*c0);
					alpha = min(1, alpha);
					alpha = max(0, alpha);
				end
			else
				alpha = 1/0.9;
				del_obj = -999;
				while (del_obj <=0) & (alpha > 0.0096)
					alpha = alpha * 0.9;
					del_obj = c2*alpha^2 + c1*alpha + c0;
				end

				if alpha <= 0.0096, alpha = 0; end
			end

			x = (1-alpha)*x + alpha*xos;
			alph(alphidx) = alpha;
			alphidx = alphidx + 1;

		end

		if subout
			xs(:,idx) = x;
			idx = idx + 1;
		end
	end

	if chat, printf('Range %g %g', min(x), max(x)), end
	if ~subout, xs(:,iter) = x; end
end

%figure, plot(alph, 'o-')

function obj = trl_obj(alpha, xinc, xos, linc, los, yi, bi, ri, R)

x = alpha*xos + (1 - alpha)*xinc;
li = alpha*los + (1 - alpha)*linc;
yp = bi .* exp(-li) + ri;
like = sum(yi .* log(yp) - yp);
penal = sum(R.penal(R, x));
obj = like - penal;
obj = -obj;
