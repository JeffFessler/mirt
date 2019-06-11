 function A = Glinear(varargin)
%function A = Glinear(options)
%|
%| Generate a geometric system matrix for tomographic projection
%| based on simple linear interpolation.
%| Contains exactly two a_ij values per pixel per projection angle.
%| This system model is very inadequate for reconstructing
%| real tomographic data.  Maybe it is handy for quick simulations.
%|
%| options
%|	nx,ny		image size
%|	nb,na		sinogram size (n_radial_bins * n_angles).
%|	ray_pix		ray_spacing / pixel_size (usually 1 or 0.5)
%|	mask [nx ny]	which pixels are to be reconstructed
%|	chat		verbosity
%| out
%|	A [nb*na nx*ny]	all of A, regardless of mask
%|			because entire size is needed for .wtf saving.
%|
%| Caller must do A = A(:,mask(:)) for masked reconstructions.
%|
%| Copyright Apr 2000, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), Glinear_test, return, end

% defaults
arg.nx = 32;
arg.ny = [];
arg.nb = [];
arg.na = [];
arg.ray_pix = 1;
arg.mask = [];
arg.chat = false;

arg = vararg_pair(arg, varargin);

if isempty(arg.ny), arg.ny = arg.nx; end
if isempty(arg.nb), arg.nb = arg.nx; end
if isempty(arg.na), arg.na = floor(arg.nb * pi/2); end
if isempty(arg.mask), arg.mask = true(arg.nx,arg.ny); end

A = Glinear_do(arg.nx, arg.ny, arg.nb, arg.na, arg.ray_pix, arg.mask, arg.chat);

function A = Glinear_do(nx, ny, nb, na, ray_pix, mask, chat);

if ray_pix < 0.5
	warning('small ray_pix will give lousy results!')
end

%
% pixel centers
%
x = [0:nx-1] - (nx-1)/2;
y = (-1)*([0:ny-1] - (ny-1)/2); % trick: to match aspire
[x,y] = ndgrid(x, y);
x = x(mask(:));
y = y(mask(:));
np = length(x); % sum(mask(:)) - total # of support pixels

angle = [0:na-1]'/na * pi;
tau = cos(angle) * x' + sin(angle) * y'; % [na,np] projected pixel center
tau = tau / ray_pix; % account for ray_spacing / pixel_size
tau = tau + (nb+1)/2; % counting from 1 (matlab)
ibl = floor(tau); % left bin
val = 1 - (tau-ibl); % weight value for left bin

ii = ibl + [0:na-1]'*nb*ones(1,np); % left sinogram index

good = ibl(:) >= 1 & ibl(:) < nb; % within FOV cases
if any(~good), warning 'FOV too small', end

%nc = np; jj = 1:np; % compact A
nc = nx * ny; jj = find(mask(:))'; % all-column A
jj = jj(ones(1,na),:);

val1 = 1-val;
if 0 % make precision match aspire?
	val = double(single(val));
	val1 = double(single(val1));
end

A1 = sparse(ii(good), jj(good), val(good), nb*na, nc); % left bin
A2 = sparse(ii(good)+1, jj(good), val1(good), nb*na, nc); % right bin
A = A1 + A2;

if 0
%	subplot(121), im(embed(sum(A)', mask)) % for compact
	subplot(121), im(reshape(sum(A), nx, ny)) % for all-column
	subplot(122), im(reshape(sum(A'), nb, na))
end

% Glinear_test()
% test demo
%
function Glinear_test
nx = 32; ny = nx; nb = 40; na = 42;
x = shepplogan(nx, ny, 1);
ix = [-(nx-1)/2:(nx-1)/2]';
iy = [-(ny-1)/2:(ny-1)/2];
rr = sqrt(outer_sum(ix.^2, iy.^2));
mask = rr < nx/2-1;
A = Glinear('nx', nx, 'ny', ny, 'nb', nb, 'na', na, 'mask', mask);
y = A * double(x(:)); % forward projection
y = reshape(y, nb, na); % reshape into sinogram array
sino = zeros(nb, na); sino(nb/2, 10) = 1;
b = reshape(A' * sino(:), nx, ny);
if im
	im plc 2 2
	im(mask, 'support mask')
	im(x, 'test image')
	im(y, 'sinogram')
	xlabel ib, ylabel ia
	im(b, 'backproject 1 ray')
end
