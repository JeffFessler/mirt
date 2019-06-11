 function y = zero_tiny_negative(x, tol)
%function y = zero_tiny_negative(x, tol)
% set negatives in y to zero, provided they are tiny negatives

if nargin < 2, tol = 100; end

xmax = max(x(:));
xmin = min(x(:));
if xmax <= 0, error 'need positive maximum', end
if xmin < -tol*eps*xmax, error 'too negative', end
y = max(x, 0);
