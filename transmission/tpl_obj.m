 function [obj, like, penal] = tpl_obj(xs, G, yi, bi, ri, R, mask, varargin)
%function [obj, like, penal] = tpl_obj(xs, G, yi, bi, ri, R, mask, [options])
%
% compute transmission penalized Poisson likelihood for each column of x
%
% model: Y_i ~ Poisson(b_i exp(-[G x]_i) + r_i)
% in
%	xs	[np,niter]	image estimates
%	G	[nd,np]		system matrix
%	yi	[nd]		transmission sinogram
%	bi	[nd]		blank scan factors
%	ri	[nd]		background (randoms, scatter, etc)
%	bi,ri:	optional (can use empty matrices)
%	yi,bi,ri must have identical dimensions
%	R			penalty object (see Robject.m)
%	mask	[np]		logical
% option
%	'use_ll' 1		do not use -(Kullback-Leibler divergence)
%				but instead use traditional log-likelihood.
%				(use_kl may be numerically inaccurate)
% out
%	obj	[niter]		objective function
%
% Copyright May 1999, Jeff Fessler

if nargin < 6, help(mfilename), error(mfilename), end
if ~isvar('bi') || isempty(bi)
	bi = ones(size(yi(:)));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi(:)));
end
if ~isvar('R'), R = []; end
if ~isvar('mask') || isempty(mask)
	mask = [];
end

arg.use_ll = 0;
arg = vararg_pair(arg, varargin);

if ~isempty(mask)
	xs = reshape(xs, size(xs,1)*size(xs,2), size(xs,3));
	xs = xs(mask(:), :);
end

like = trl_like(xs, G, yi, bi, ri, arg.use_ll); % likelihood

penal = zeros(size(like));
if isempty(R)
	warning 'no penalty'

% trick: if the "R" input is sparse, then it is actually the "C" matrix
elseif issparse(R)
	for kk=1:ncol(xs)
		penal(kk,1) = sum((R * xs(:,kk)).^2)/2;
	end

elseif isstruct(R)
	for kk=1:ncol(xs)
		penal(kk,1) = R.penal(R, xs(:,kk));
	end
else
	error 'R'
end

obj = like - penal;

%printf('%g\t', [like penal obj]')


% Poisson transmission likelihood for each column of x
function like = trl_like(xs, G, yi, bi, ri, use_ll)

like = zeros(ncol(xs),1);
for kk=1:ncol(xs)
	x = xs(:,kk);
	li = G * x;
	yp = bi .* exp(-li) + ri;	% predicted measurements
	if use_ll
		if any((yp == 0) & (yi > 0)), error 'model mismatch', end
		like(kk) = sum(yi .* log(yp) - yp);
	else
		yid = yi; yid(yi == 0) = 1;
		ypn = yp; ypn(yi == 0) = 1;
		like(kk) = sum(yi .* log(ypn ./ yid) + (yi - yp));
	end
end
