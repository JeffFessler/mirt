% ct_beta_study
% study how to adjust beta when pixel size changes in X-ray CT recon.
% results:
% changing nx has neglible effect, as expected
% changing pixel size has minimal effect when fwhm is converted to mm!
% values used CT studies:
% l2b=9, fov=500, fwhm = 1.62 mm
% l2b=9, fov=250, fwhm = 1.63 mm
% so we can  use the same beta regardless of pixel size and get nearly same fwhm
% caution: if other factors change, such as the # of views, may need to revisit!

%nx = 256;
%nx = 128;
nx = 64;
ny = nx;

if ~isvar('R')
	mask = true(nx,ny);
	l2b = 9;
	R = Robject(mask, 'edge_type', 'tight', 'beta', 2^l2b);
end

pixel_size_list = [500/512 250/512 0.25];
for ii=1:length(pixel_size_list)
	pixel_size = pixel_size_list(ii);

	sys = ct_sys('nx', nx, 'ny', ny, 'pixel_size', pixel_size);

	G = Gtomo2_dscmex(sys, 'nthread', 2);

%	psf = qpwls_psf(G, C, 2^l2b, mask, 1);
	psf = qpwls_psf(G, R, 1, mask, 1);
	fwhm = fwhm2(psf);
	printf('fwhm in mm = %g', fwhm * pixel_size)
end
