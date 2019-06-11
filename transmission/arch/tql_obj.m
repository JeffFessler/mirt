 function obj = tql_obj(xs, G, yi, bi, ri, C, mask)
%function obj = tql_obj(xs, G, yi, bi, ri, C, mask)
%	compute quadratically penalized Poisson likelihood for each column of x
%
%	Copyright May 1999, Jeff Fessler

if nargin < 3, help(mfilename), error(mfilename), end
if (nargin < 4 || isempty(bi))
	bi = ones(size(yi(:)));
end
if (nargin < 5 || isempty(ri))
	ri = zeros(size(yi(:)));
end
if (nargin < 6 || isempty(C))
	C = 0;
end
if nargin == 7
	xs = reshape(xs, size(xs,1)*size(xs,2), size(xs,3));
	xs = xs(find(mask(:)), :);
end

obj = zeros(ncol(xs),1);
for kk=1:ncol(xs)
	x = xs(:,kk);
	yp = bi .* exp(-(G * x)) + ri;	% predicted measurements
	obj(kk) = sum(yi .* log(yp) - yp) - norm(C * x).^2 / 2;
end
