% spect_1d_spot_example.m
%
% Investigate the "hole in sphere" effect using a 1D model (for speed).
% This simulates 1D rect functions of various widths on a uniform background,
% seen through a wide 1D blur.  When we run ML-EM, we find that the ML estimate
% exhibits over shoot or undershoot at the center of the ROI, depending on
% the size of the ROI relative to the width of the reconstructed PSF.
%
% Copyright 2006-5-31, Jeff Fessler, The University of Michigan

% true emission image: collection of 1D rect functions of various widths
% on a uniform background
if ~isvar('xtrue'), printm 'xtrue'
	over = 2^4;
	ig = image_geom('nx', over*256, 'ny', 100, 'dx', 1/over);
	[xx yy] = ndgrid(ig.x, [1:ig.ny]/ig.ny);
	fwhm = 70*yy;
	tmp = double6(abs(xx - 0.15*ig.fov) < fwhm/2);
	fwhm = fwhm(1,:)';
	xtrue = abs(xx) < ig.fov * 0.4;
	xtrue = xtrue + 2 * tmp;
	xtrue = downsample2(xtrue, [over 1]);

	ig = image_geom('nx', ig.nx/over, 'ny', ig.ny, 'dx', 1);
	xtrue(round(ig.nx*0.3),:) = 3;

	ixpeak = imax(xtrue(:,1) .* (ig.x > 0));
	im plc 2 3
	im subplot 4
	if im, plot(1:ig.ny, fwhm1(xtrue, 'imid', ixpeak)), end
	ax = {ig.x, fwhm};
	if im, im(1, ax{:}, xtrue, 'xtrue'), axis normal, cbar, end
prompt
end

% system model: 1D blur
if ~isvar('G'), printm 'G'
	psf = gaussian_kernel(12);
	if im
		im subplot 5, plot(psf, '-o')
		title(sprintf('1D System PSF.  FWHM=%g', fwhm1(psf)))
	end
	G = Gblur(ig.mask, 'psf', psf);
prompt
end

% simulate noiseless blurry data
if ~isvar('yt'), printm 'yt'
	yt = G * xtrue;
	if im, im(2, ax{:}, yt, 'yt'), axis normal, cbar, end
prompt
end

% ordinary EM
if ~isvar('xmlem'), printm 'xmlem'
	f.niter = 60+1;
	xinit = ig.ones;
	xmlem = eml_em(xinit(ig.mask), G, yt(:), 1, 0, [], f.niter);
	xmlem = xmlem(:,end);
	xmlem = ig.embed(xmlem);
	if im, im(3, ax{:}, xmlem, 'x ML-EM'), axis normal, cbar, end
prompt
end

if im
	im subplot 4
	ix = imax(xmlem(:,1) .* (ig.x < 0));
	ix = ix + [-30:30];
	plot(ix, xmlem(ix,1)-1)
	axis tight
	xlabel 'pixel', title 'local psf'
end

if im
	im subplot 6
	plot( ...
		fwhm, max(xmlem, [], 1), 'y-', ...
		fwhm, xmlem(ixpeak,:), 'c-' ...
	)
	%	fwhm, min(xmlem, [], 1), 'g-', ...
	ir_legend({'roi max', 'roi center'})
	axis tight
	xlabel 'width'
	grid

	im subplot 5
	ii = imin(abs(fwhm - 29));
	plot(	ig.x, xmlem(:,ii), 'c-', ...
		ig.x, xtrue(:,ii), 'y-' ...
		)
	axis tight
	title(sprintf('"profile" for object width=%g', fwhm(ii)))
end
