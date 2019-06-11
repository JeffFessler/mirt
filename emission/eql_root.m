 function x = eql_root(a, b, c)
%function x = eql_root(a, b, c)
% Numerically stable method for computing the positive root
% of the quadratic polynomial: -ax^2 -2bx + c, with a >= 0.
% This polynomial arises in quadratically-penalized likelihood
% methods for the Poisson emission image reconstruction problem.
%
% Copyright Apr. 2000, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

if any(a < 0), error 'need a > 0', end

x = zeros(size(a));

j = a == 0;
x(j) = c(j) ./ b(j) / 2;		% 1st-order case

det = sqrt(b.^2 + a .* c);		% determinant / 2

j = (a > 0) & (b > 0);
x(j) = c(j) ./ (det(j) + b(j));

j = (a > 0) & (b <= 0);
x(j) = (det(j) - b(j)) ./ a(j);
