 function cost = qpwls_cost(xs, G, W, yy, C, mask)
%function cost = qpwls_cost(xs, G, W, yy, C, mask)
% compute QPWLS cost for each column of x
%
% Copyright Apr 1999, Jeff Fessler

if nargin < 3, help(mfilename), error(mfilename), end

if nargin == 6
	xs = reshape(xs, size(xs,1)*size(xs,2), size(xs,3));
	xs = xs(find(mask(:)), :);
end

cost = zeros(ncol(xs),1);
for kk=1:ncol(xs)
	x = xs(:,kk);
	resid = yy - G * x;	% predicted measurements
	cost(kk) = resid' * W * resid / 2 + norm(C * x).^2 / 2;
end

cost = reale(cost);	% trick: x'*x is not real for complex values
