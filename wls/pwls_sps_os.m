 function [xs, info] = pwls_sps_os(x, yi, wi, Ab, R, niter, pixmax, ...
		denom, aai, relax0, chat)
%function [xs, info] = pwls_sps_os(x, yi, wi, Ab, R, niter, pixmax, ...
%|		denom, aai, relax0, chat)
%|
%| penalized weighted least squares estimation/reconstruction
%| using separable paraboloidal surrogates algorithm with
%| (optionally relaxed) ordered subsets.  (relaxation ensures convergence.)
%|
%| cost(x) = (y-Gx)' W (y-Gx) / 2 + R(x)
%|
%| in
%|	x	[np 1]		initial estimate
%|	yi	[nb na]		measurements (noisy sinogram)
%|	wi	[nb na]		weighting sinogram (or [] for uniform)
%|	Ab	[nd np]		Gblock object, aij >= 0 required!
%|					or sparse matrix (implies nsubset=1)
%|	R			penalty object (see Reg1.m)
%|	niter			# of iterations (including 0)
%|
%| optional
%|	pixmax	[1] or [2]	max pixel value, or [min max] (default [0 inf])
%|	denom	[np 1]		precomputed denominator
%|	aai	[nb na]		precomputed row sums of |Ab|
%|	relax0	[1] or [2]	relax0 or (relax0, relax_rate)
%|	chat
%|
%| out
%|	xs	[np niter]	iterates
%|	info	[niter 1]	time
%|
%| Copyright 2002-2-12, Jeff Fessler, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

cpu etic
info = zeros(niter,1);

Ab = block_op(Ab, 'ensure'); % make it a block object (if not already)
nblock = block_op(Ab, 'n');
starts = subset_start(nblock);

if ~isvar('niter')	|| isempty(niter),	niter = 1;	end
if ~isvar('pixmax')	|| isempty(pixmax),	pixmax = inf;	end
if ~isvar('chat')	|| isempty(chat),	chat = false;	end
if isempty(wi)
	wi = ones(size(yi));
end
if ~isvar('aai') || isempty(aai)
	aai = reshape(sum(Ab'), size(yi)); % a_i = sum_j |a_ij|
					% requires real a_ij and a_ij >= 0
end

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

if length(pixmax) == 2
	pixmin = pixmax(1);
	pixmax = pixmax(2);
elseif length(pixmax) == 1
	pixmin = 0;
else
	error pixmax
end

%
% likelihood denom, if not provided
%
if ~isvar('denom') || isempty(denom)
	denom = Ab' * col(aai .* wi);	% requires real a_ij and a_ij >= 0
end, clear aai

if ~isvar('R') || isempty(R)
	pgrad = 0;		% unregularized default
	Rdenom = 0;
end

[nb na] = size(yi);


%
% loop over iterations
%

xs = zeros(numel(x), niter, class(x));
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

		li = Ab{iblock} * x;
		li = reshape(li, nb, length(ia));
		resid = wi(:,ia) .* (yi(:,ia) - li);
		grad = Ab{iblock}' * resid(:); % G' * W * (y - G*x)

		if ~isempty(R)
			pgrad = R.cgrad(R, x);
			Rdenom = R.denom(R, x);
		end

		num = nblock * grad - pgrad;
		den = denom + Rdenom;

		x = x + relax * num ./ den;	% relaxed update
		x = max(x,pixmin);		% lower bound
		x = min(x,pixmax);		% upper bound
	end

	if chat, printm('range(x) = %g %g', min(x), max(x)), end
	xs(:,iter) = x;
	info(iter,1) = cpu('etoc');
end
