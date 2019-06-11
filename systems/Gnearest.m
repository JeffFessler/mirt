 function G = Gnearest(nx, ny, nb, na, mask)
%|function G = Gnearest(nx, ny, nb, na, mask)
%|
%| Generate a geometric system matrix for tomographic projection
%| based on simple nearest-neighbor "interpolation."
%| This system model is completely inadequate for reconstructing
%| real tomographic data, but is marginally useful for simple simulations.
%|
%| The image size is nx * ny.
%| The sinogram size is nb * na (n_radial_bins X n_angles).
%| Returned G is [nb * na, nx * ny] - this size is need for .wtf saving.
%|
%| Jeff Fessler

if nargin < 1, nx = 64; end
if nargin < 2, ny = nx; end
if nargin < 3, nb = nx; end
if nargin < 4, na = floor(nb * pi/2); end
if nargin < 5, mask = true(nx,ny); end

warning 'this function is obsolete.  use Gtomo2_strip() instead'

% pixel centers
[x,y] = ndgrid([0:nx-1] - (nx-1)/2, [0:ny-1] - (ny-1)/2);
x = x(mask(:));
y = y(mask(:));
np = length(x);		% sum(mask(:)) - total # of support pixels

angle = [0:na-1]'/na * pi;
tau = cos(angle) * x' + sin(angle) * y';	% [na,np] projected pixel center
ii = round(tau + (nb-1)/2);			% counting from 1 (matlab)
good = ii(:) >= 1 & ii(:) <= nb;
if any(~good), warning 'FOV too small', end

ii = ii + [0:na-1]' * nb * ones(1,np);

%np = sum(mask(:));
%nc = np;	jj = 1:np;		% compact G
nc = nx * ny;	jj = find(mask(:))';	% all-column G
jj = jj(ones(1,na),:);

G = sparse(ii(good), jj(good), ones(sum(good),1), nb*na, nc);

if 0
%	subplot(121), im(embed(sum(G)', mask))		% for compact
	subplot(121), im(reshape(sum(G), nx, ny))	% for all-column
	subplot(122), im(reshape(sum(G'), nb, na))
end
