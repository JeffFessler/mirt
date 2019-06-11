% psfs_uhe(file_psfs, radius, dim, ny, sxz, sy)
%
%   Generate ultra-high-energy collimator PSF's for SPECT 3D reconstruction
%   where 
%   file_psfs:	Output file the PSF's are to be written to ([] for no file)
%   radius:	Object-to-detector distance (mm)
%   dim:	Dimension of generated PSF's, either 1 or 2 (default is 2)
%   ny: 	Number of slices in y direction (default is 64)
%   sxz:	Pixel size (mm) in directions x-z parallel to detector (default is 7.196)
%   sy:		Pixel size (mm) in direction y perpendicular to detector (default is equal to sxz)
%
% Qiang (Victor) Lin
% Modified by A. Yendiki, 2/5/02

function psfs = psfs_uhe(file_psfs, radius, dim, ny, sxz, sy)

if ~isvar('dim') | isempty(dim)
	dim = 2;	% Default is 2-D PSF's
end
if ~isvar('ny') | isempty(ny)
	ny = 64;	% Default slice # in y direction
end
if ~isvar('sxz') | isempty(sxz)
	sxz = 7.196;	% Default pixel size in x and z (mm)
end
if ~isvar('sy') | isempty(sy)
	sy = sxz;	% Default pixel size in y (mm)
end

% Slice-to-detector distance (mm)
distance = radius + [-(ny-1)/2:(ny-1)/2] * sy;
fwhm = getFWHM_uhe(distance);

% Determine the size of the final matrix (pixels)
nx_pad = 2 * ceil(fwhm(end) / sxz) + 1;
nz_pad = nx_pad;

% Initialize the matrix
if dim == 2
	psfs = zeros(nx_pad, nz_pad, ny);
elseif dim == 1
	psfs = zeros(ceil(nx_pad/2), ny);
else
	error('Dimension of generated PSF''s can only be 1 or 2')
end

% Points where Gaussian will be sampled
x = [-(nx_pad-1)/2 : (nx_pad-1)/2]' * sxz;

% Generate PSF at each distance from the detector
for iy=1:ny
	% Slice-to-detector distance (mm)
	if distance(iy) >= 0
		% Standard deviation of Gaussian (mm)
		sigma = fwhm(iy) / sqrt(log(256));

		% Generate 1D Gaussian
		tmp = normcdf(x+sxz/2, 0, sigma) - normcdf(x-sxz/2, 0, sigma);

		if dim == 2
			% Generate 2D Gaussian (separable)
			if nz_pad > 1
				tmp = tmp * tmp';
			end
		end
		tmp = tmp / sum(tmp(:));       % Scale to get unity sum
		tmp(find(tmp < 0.000001)) = 0.0;

		% Save PSF in corresponding slice of 3-D PSF matrix
		if dim == 2
			psfs(:,:,iy) = tmp;
		else
			% Keep right part of 1D Gaussian
			psfs(:,iy) = tmp(ceil(nx_pad/2) : end);
		end
		%disp(psfs(:,:,iy))
		%pause
	end
end
%disp(tmp)       % See how close to unity the last one is
%im(psfs)

if ~isempty(file_psfs)
	fld_write(file_psfs, psfs)
end

