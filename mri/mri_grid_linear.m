  function [xhat, yhat, xg, kg] = mri_grid_linear(kspace, ydata, N, fov)
%|function [xhat, yhat, xg, kg] = mri_grid_linear(kspace, ydata, N, fov)
%| very crude "gridding" based on linear interpolation.
%| not recommended as anything but a straw man or perhaps
%| for initializing iterative methods.
%| in
%|	kspace		[M 2]	kx and ky k-space sample locations
%|	ydata		[M 1]	complex Fourier sample values
%|	N			desired image size
%|	fov			field of view
%|		kspace and fov must be in reciprocal units!
%| out
%|	xhat		[N N]	image estimate
%|	yhat		[N N]	gridding k-space data
%|	xg	{[N1 1] [N2 1]}	object-space uniform sample locations
%|	kg	{[N1 1] [N2 1]}	k-space uniform sample locations in 1d
%|
%| Copyright 2003-7-23, Jeff Fessler, University of Michigan

if nargin == 1 && streq(kspace, 'test'), mri_grid_linear_test, return, end
if nargin < 4, help(mfilename), error(' '), end

if length(N) == 1, N = [N N];, end
if length(N) ~= 2, error 'bad N', end
if length(fov) == 1, fov = [fov fov];, end
if length(fov) ~= 2, error 'bad fov', end

kg{1} = [-N(1)/2:N(1)/2-1]/fov(1);
kg{2} = [-N(2)/2:N(2)/2-1]/fov(2);
[k1gg, k2gg] = ndgrid(kg{1}, kg{2});
yhat = griddata(kspace(:,1), kspace(:,2), ydata, k1gg, k2gg, 'linear');
yhat(isnan(yhat)) = 0;

xg{1} = [-N(1)/2:N(1)/2-1]'/N(1) * fov(1);
xg{2} = [-N(2)/2:N(2)/2-1]'/N(2) * fov(2);
dk = 1 ./ fov;
xhat = fftshift(ifft2(fftshift(yhat))) * prod(dk) * prod(N);
tmp1 = nufft_sinc(xg{1} / fov(1));
tmp2 = nufft_sinc(xg{2} / fov(2));
xhat = xhat ./ (tmp1 * tmp2'); % gridding post-correction for linear interp


%
% test function
%
function mri_grid_linear_test

fov = 256; % [mm] typical brain FOV
N0 = 64; % nominal image size

t = linspace(0, N0/2*2*pi, N0^2)'; % crude spiral:
kspace = N0/2*(1/fov)*[cos(t) sin(t)] .* (t(:,[1 1]) / max(t));

Ndisp = 256; % display images like this...
x1d = [-Ndisp/2:Ndisp/2-1] / Ndisp * fov;
[x1dd x2dd] = ndgrid(x1d, x1d);

obj = mri_objects('case1');
xtrue = obj.image(x1dd, x2dd);
ytrue = obj.kspace(kspace(:,1), kspace(:,2));

im plc 2 3
clim = [0 2];
im(1, x1d, x1d, xtrue, 'x true', clim), cbar

Ng = 128;
[xhat yhat xg kg] = mri_grid_linear(kspace, ytrue, Ng, fov);

im(2, kg{1}, kg{2}, abs(yhat), '|y_{grid}|'), cbar
im(3, xg{1}, xg{2}, abs(xhat), '|x| "gridding"', clim), cbar
