% em_wls_test_setup.m
%
% create sample image, system matrix, and sinograms for examples
% and testing of WLS algorithm for emission tomography.
% creates: sg ig xtrue mask G ci yi wi
%
% Copyright June 2000, Jeff Fessler, University of Michigan

if ~isvar('ig')
	ig = image_geom('nx', 64, 'ny', 60, 'fov', 500);
	ig.mask = ig.circ > 0; % support mask (which pixels do we estimate?)
end
if ~isvar('sg')
	sg = sino_geom('par', 'nb', ig.nx+2, 'dr', 528 / (ig.nx+2), ...
		'strip_width', 'dr');
end


% true emission image and attenuation map
if ~isvar('xtrue')
	xtrue = read_zubal_emis('nx', ig.nx, 'ny', ig.ny);
	mumap = read_zubal_attn('nx', ig.nx, 'ny', ig.ny);
	im pl 3 3
	im(1, xtrue, 'true emission image'), cbar
	im(2, mumap, 'attenuation map'), cbar
	im(3, ig.mask, 'support mask'), cbar
end


% system matrix G
if ~isvar('G'), printm 'G'
	G = Gtomo2_strip(sg, ig, 'single', 1);

	if isvar('f.wtf') && has_mex_jf
		if exist(f.wtf), delete(f.wtf), end
		wtf_write(f.wtf, G, ig.nx, ig.ny, sg.nb, sg.na);
	end
	if isvar('f.wtr') && has_aspire && has_mex_jf
		if exist(f.wtr), delete(f.wtr), end
		os_run(sprintf('wt -chat 0 col2row %s %s', f.wtr, f.wtf))
	end
%	printm('condest(G'*G) = %g', condest(G'*G))
end


% for default 64 size, l2b=9 yields fwhm = 1.46 pixel or 11.4 mm
if 0, printm 'explore FWHM of PSF'
	f.l2b = 9;
	Rq = Reg1(ig.mask, 'beta', 2^f.l2b);
	psf = qpwls_psf(G, Rq.C, 1, ig.mask, 1, 'dx', ig.dx);
	im(7, psf, 'PSF')
return
end


% noisy Poisson sinogram measurements
if ~isvar('yi'), printm 'yi'
	proj = G * xtrue;
	ci = G * mumap;
	printm('attenuation length range %g %g', min(ci(:)), max(ci(:)))
	ci = exp(-ci);
	if isvar('f.count')
		count = f.count;
	else
		count = 1e5;
	end
	ci = count / sum(col(ci .* proj)) * ci;
	ci = dsingle(ci);
	ytrue = ci .* proj;
	if ~isvar('randpercent')
		randpercent = 10;
	end
	ri = randpercent / 100 * mean(ytrue(:)) * ones(size(ytrue));
	ri = dsingle(ri);
	rng(0)
	ypi = poisson(ytrue + ri);

	im(4, ytrue, 'ytrue: true projections'), cbar
	im(5, ypi, 'ypi: noisy projections'), cbar
	clear count randpercent ytrue proj

	% precorrect measurements and compute weights
	yi = (ypi - ri) ./ ci;
	wi = ci.^2 ./ max(ypi, 1);
	clear ypi ri ci
end


% FBP reconstruction
if ~isvar('f.fbp_window')
	f.fbp_window = [];
end
if ~isvar('xfbp'), printm 'fbp'
	tmp = fbp2(sg, ig, 'window', f.fbp_window);
	xfbp = fbp2(max(yi,0), tmp);
	im(6, max(xfbp,0), 'FBP Reconstruction'), cbar
	clear tmp
end

if isvar('f.yi')
	fld_write(f.yi, yi, 'check', 0)
end
if isvar('f.wi')
	fld_write(f.wi, wi, 'check', 0)
end
if isvar('f.mask')
	fld_write(f.mask, mask, 'check', 0)
end
