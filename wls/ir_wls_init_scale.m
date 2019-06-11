 function [xnew, scale] = ir_wls_init_scale(A, y, xold)
%function [xnew, scale] = ir_wls_init_scale(A, y, xold)
%|
%| Find the scaled version of the image "x" that best fits the data.
%| Despite the name, it currently is "LS" not "WLS"
%|
%| out
%|	xnew = scale * xold
%|	scale = argmin_s |y - s A x|_2^2
%|
%| 2014, Jeff Fessler, University of Michigan

if nargin == 1 && streq(A, 'test'), ir_wls_init_scale_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

if nargin < 3
	xold = A' * y; % default that is sensible only for some applications
end
tmp = A * xold;
scale = sum(col(conj(tmp) .* y), 'double') / sum(col(abs(tmp).^2), 'double');
xnew = scale * xold;


% ir_wls_init_scale_test
% test with complex data!
function ir_wls_init_scale_test
x = [1:5]' + 1i * [2:6]';
rng(0)
A = randn(8,5) + 1i * randn(8,5);
scale1 = 3;
y = scale1 * A * x;
[x2 scale2] = ir_wls_init_scale(A, y, x);
equivs(scale1, scale2)
