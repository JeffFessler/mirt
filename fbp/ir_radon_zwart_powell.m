 function output = ir_radon_zwart_powell(theta, rr)
%function output = ir_radon_zwart_powell(theta, yr)
%|
%| Compute analytic 2D Radon transform of Zwart-Powell box spline.
%|
%| in
%| theta	ray angle in radian
%| rr	distance between the point and the ray (normalized by the pixel size)
%| output: radon transform of the Zwart-Powell box spline element
%|
%| This is the modified version of the code written by [1]
%| to avoid symbolic math operation.
%| Code written by Seongjin Yoon, Univ. of Michigan, Jan 2015
%|
%| Reference
%| [1] A. Entezari, M. Nilchian, and M. Unser, "A box spline calculus
%| for the discretization of computed tomography reconstruction problems,"
%| IEEE Trans. Med. Imaging, vol. 31, no. 8, pp. 1532–1541, 2012.
%|
%| 2015-08-10 Jeff Fessler, added self test and parallelized

if nargin == 1 && streq(theta, 'test'), ir_radon_zwart_powell_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

dim = size(theta);
theta = theta(:);

zeta = zeros(numel(theta), 4);
zeta(:,1) = cos(theta);
zeta(:,2) = sin(theta);
zeta(:,3) = zeta(:,1) + zeta(:,2);
zeta(:,4) = zeta(:,1) - zeta(:,2);

cond = abs(zeta) >= eps('single');
N = sum(cond,2);

output = BoxSp4(rr(:), zeta, cond, N) ./ factorial(N-1);
output = reshape(output, dim);


function output = BoxSp0(y, N)
%output = heaviside(y) .* y.^(N-1);
if any(size(y) ~= size(N)), keyboard, end
output = (y >= 0) .* y.^(N-1);


function output = BoxSp1(y, zeta, cond, N)
good = cond(:,1);
output = (BoxSp0(y+0.5*zeta(:,1), N) ...
	- BoxSp0(y-0.5*zeta(:,1), N)) ./ zeta(:,1);
output(~good) = BoxSp0(y(~good), N(~good));


function output = BoxSp2(y, zeta, cond, N)
good = cond(:,2);
output = (BoxSp1(y+0.5*zeta(:,2), zeta, cond, N) ...
	- BoxSp1(y-0.5*zeta(:,2), zeta, cond, N)) ./ zeta(:,2);
output(~good) = BoxSp1(y(~good), zeta(~good,:), cond(~good,:), N(~good));


function output = BoxSp3(y, zeta, cond, N)
good = cond(:,3);
output = (BoxSp2(y+0.5*zeta(:,3), zeta, cond, N) ...
	- BoxSp2(y-0.5*zeta(:,3), zeta, cond, N)) ./ zeta(:,3);
output(~good) = BoxSp2(y(~good), zeta(~good,:), cond(~good,:), N(~good));


function output = BoxSp4(y, zeta, cond, N)
good = cond(:,4);
output = (BoxSp3(y+0.5*zeta(:,4), zeta, cond, N) ...
	- BoxSp3(y-0.5*zeta(:,4), zeta, cond, N)) ./ zeta(:,4);
output(~good) = BoxSp3(y(~good), zeta(~good,:), cond(~good,:), N(~good));


% ir_radon_zwart_powell_test()
function ir_radon_zwart_powell_test
theta = linspace(0,pi,181);
r = linspace(-1,1,101) * 2;
[tt rr] = ndgrid(theta, r);
sino = ir_radon_zwart_powell(tt, rr);
im('colorneg', r, theta, sino')
