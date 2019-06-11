 function [xs, ni] = eml_sps(x, Gt, yi, ci, ri, niter, curv)
%function [xs, ni] = eml_sps(x, Gt, yi, ci, ri, niter, curv)
%	One iteration of the ML-SPS algorithm for emission Poisson problem
%	(separable paraboloidal surrogates)
%	model: Y_i ~ Poisson(c_i [G x]_i + r_i)		WITH r_i > 0 REQUIRED!
%	in:
%		Gt	transpose of system matrix
%		see em_fbp.m for model, G, yi, ci, ri
%	out:
%		x [np,niter]	updated image vectors each iteration
%
%	Copyright Mar 2000, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

[nb, na] = size(yi);

if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end
if ~isvar('niter') || isempty(niter)
	niter = 1;
end

eml_check(yi, ci, ri)

	gi = sum(Gt)';	% g_i = sum_j g_ij

xs = zeros(numel(x), niter);
xs(:,1) = x(:);
%
%	loop over iterations
%
for ii=2:niter
	li = reshape(Gt' * x(:), size(yi));	% l=G*x "line integrals"
	yb = ci .* li + ri;			% predicted measurement means 

	%	curvatures
	ni = eml_curvature(yi, ci, ri, li, yb, curv);

	dothi = ci .* (yi ./ yb - 1);  

	x = x + (Gt * dothi(:)) ./ (Gt * (gi .* ni(:)));
	x = max(x,0);

	xs(:,ii) = x;
end
