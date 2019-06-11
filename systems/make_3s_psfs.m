 function psfs = make_3s_psfs(ny, sx, obj2det, fwhm0, slope)
%function psfs = make_3s_psfs(ny, sx, obj2det, fwhm0, slope)
% generate Gaussian PSFs for SPECT 3D testing/simulation
% in
%	ny		# of image pixels in one direction
%	sx		# pixel size [mm]
%	obj2det		object-to-detector distance [mm]
%	fwhm0		FHWM(0) [mm]
%	slope		slope of FHWM(z) (unitless)
% out
%	psfs		[npsf npsf ny]

if nargin < 5
	help(mfilename), error(mfilename)
end

%
% system-dependent factors
%
if 0
	ny = 64;
	sx = 0.72;		% pixel size in cm
	obj2det = 26;		% detector orbit radius of rotation
	fwhm0 = 0.85611;	% collimator FWHM at z=0
	fwhm30 = 2.83613;
	slope = sqrt(fwhm30^2 - fwhm0^2) / 30;
end

y = -( [0:(ny-1)]-(ny-1)/2 ) * sx;
z = obj2det - y;	% distance of each plane from detector
fwhm = sqrt(fwhm0^2 + (z * slope).^2);

nx_pad = 2 * ceil(max(fwhm/sx)) + 1;
nz_pad = nx_pad;	% use this for 2D blur
%nz_pad = 1;		% use this for 1D blur

psfs = zeros(nx_pad, nz_pad, ny);

pn = jf_protected_names;

x = [-(nx_pad-1)/2 : (nx_pad-1)/2]' * sx;
for iy=1:ny
	sig = fwhm(iy) / sqrt(log(256));
	tmp = pn.normcdf(x+sx/2, 0, sig) - pn.normcdf(x-sx/2, 0, sig);
	if nz_pad > 1
		tmp = tmp * tmp';	% 2D gaussian is separable
	else
		tmp = tmp;		% 1D gaussian
	end
	psfs(:,:,iy) = tmp / sum(tmp(:));	% sum to unity
end

printf('the last PSF sums to this value: %g', sum(tmp(:)))

izc = imin(abs(z - obj2det))
if im
	clf, subplot(121), plot(z, fwhm, 'c-o', z(izc), fwhm(izc), 'y*'), grid
	xlabel 'z [cm]', ylabel 'FWHM(z)'
	im(122, psfs, 'PSFs')
end

%if strcmp(input('do you want to save these psfs? ', 's'), 'y')
%	save psfs psfs
%end
