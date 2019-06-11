% psf_mismatch_example2.m
%
% 2D example showing the effect of PSF mismatch on ML-EM algorithm
%
% Copyright 2001-8-24, Jeff Fessler, The University of Michigan

if ~has_aspire, return, end

%
% generate data
%
if ~isvar('yi'), printm 'setup psf_mismatch_example2'
	f.dir = test_dir;
	f.dsc0 = [f.dir 't0.dsc'];
	f.dsc1 = [f.dir 't1.dsc'];
	f.dsc2 = [f.dir 't2.dsc'];
	f.wtf0 = strrep(f.dsc0, 'dsc', 'wtf');
	f.wtf1 = strrep(f.dsc1, 'dsc', 'wtf');
	f.wtf2 = strrep(f.dsc2, 'dsc', 'wtf');
	if 1
		Fwhm0 = 5;
		os_run(['wt -chat 0 dsc 12 fwhm_detector 5 >! ' f.dsc0])
		os_run(['echo y | wt -chat 0 gen ' f.dsc0])
		Fwhm1 = 7;
		os_run(['wt -chat 0 dsc 12 fwhm_detector 7 >! ' f.dsc1])
		os_run(['echo y | wt -chat 0 gen ' f.dsc1])
		Fwhm2 = 2;
		os_run(['wt -chat 0 dsc 12 fwhm_detector 2 >! ' f.dsc2])
		os_run(['echo y | wt -chat 0 gen ' f.dsc2])
	end

%	G0 = Gtomo2_wtmex('f.wtf0'); % cannot use with multiple .wtf's

	[t nx ny nb na] = wtf_read(f.wtf0);
	ig = image_geom('nx', nx, 'ny', ny, 'dx', 1);
	ig.mask = reshape(sum(t) ~= 0, nx, ny);

	G0 = Gsparse(f.wtf0, 'mask', ig.mask);
	G1 = Gsparse(f.wtf1, 'mask', ig.mask);
	G2 = Gsparse(f.wtf2, 'mask', ig.mask);

	sg = sino_geom('par', 'nb', nb, 'na', na, 'dr', 1);
	xtrue = ellipse_im(ig, [6.5 -0.5 5 5 0 100], 'oversample', 4);

	ri = 1;
	yi = sg.shape(G0 * xtrue(ig.mask)) + ri;
	im(yi, 'yi'), cbar
prompt
end

if ~isvar('Gb2'), printm 'Gbs'
	f.nblock = 6;
	Gb0 = Gblock(G0, f.nblock);
	Gb1 = Gblock(G1, f.nblock);
	Gb2 = Gblock(G2, f.nblock);
prompt
end


% uniform initial image
xinit = ones(sum(ig.mask(:)),1);

% FBP
if 0 && ~isvar('xfbp'), printm 'fbp'
	xfbp = em_fbp(sg, ig, yi, 1, ri);
	im(xfbp), cbar
prompt
end

if ~isvar('x2') && 0
	x0 = eml_em(xinit, G0, yi(:), 1, ri, 'isave', 'all', 'niter', f.niter);
	x1 = eml_em(xinit, G1, yi(:), 1, ri, 'isave', 'all', 'niter', f.niter);
	x2 = eml_em(xinit, G2, yi(:), 1, ri, 'isave', 'all', 'niter', f.niter);
	x0 = ig.embed(x0);
	x1 = ig.embed(x1);
	x2 = ig.embed(x2);
end

%
% OS-EM iterations
%
if ~isvar('xo2'), printm 'run os-em'
	f.niter = 40;
	rri = ri * ones(nb,na);
	osem = @(G) eml_osem(xinit, G, ...
		yi, [], rri, 'niter', f.niter, 'isave', 'all');
	xo0 = osem(Gb0);
	xo1 = osem(Gb1);
	xo2 = osem(Gb2);
	xo0 = ig.embed(xo0);
	xo1 = ig.embed(xo1);
	xo2 = ig.embed(xo2);
	im(xo2)
prompt
end

if ~isvar('fw2'), printm 'find widths'
	ix = 39+[-9:9];
	iy = 33+[-9:9];

	xx0 = xo0(ix,iy,2:end);
	xx1 = xo1(ix,iy,2:end);
	xx2 = xo2(ix,iy,2:end);
	im(xx0)
	clear fw0 fw1 fw2
	for ii=1:f.niter
		fw0(ii) = fwhm2(xx0(:,:,ii));
		fw1(ii) = fwhm2(xx1(:,:,ii));
		fw2(ii) = fwhm2(xx2(:,:,ii));
	end

	xx = xtrue(ix,iy);
	fx = fwhm2(xx);
end

if 1 && im
	im clf
	ii = 1:f.niter;
	plot(ii(1:4:end), fw0(1:4:end), 'yo', ...
		ii(1:4:end), fw1(1:4:end), 'cx', ...
		ii(1:4:end), fw2(1:4:end), 'g+', ...
		ii, fx + 0*ii, 'b-', ...
		ii, fw0, 'y-', ...
		ii, fw1, 'c-', ...
		ii, fw2, 'g-')
	axisy(3,11)
	legend(	sprintf('True system PSF, FWHM=%g', Fwhm0), ...
		sprintf('Model PSF too big, FWHM=%g', Fwhm1), ...
		sprintf('Model PSF too small, FWHM=%g', Fwhm2), ...
		'ideal')
	xlabel 'Iteration of OS-EM algorithm'
	ylabel 'Width of reconstructed circle'
end
