% recon_limited_angle1.m
% illustrate limited angle reconstruction using edge-preserving regularization.

if ~isvar('sino')
	ig = image_geom('nx', 100, 'ny', 128, 'dx', 2);
	xell = [0 0 85 115 0 100;
		0 -60 30 20 0 -80; % mouth
		30 20 25 35 30 20; % eyes
		-30 20 25 35 -30 20;
		35 25 7 7 0 -100; % pupils
		-15 25 7 7 0 -100;
		0 75 60 15 0 -50; % hat
		]; % it is not north park, it is ...
	[xtrue xell] = ellipse_im(ig, xell, 'oversample', 3);

	im plc 2 3, clim = [0 130];
	im(1, xtrue, clim, 'xtrue'), cbar

	% noiseless sinogram for illustration
%	f.orbit = 140; % some missing views
	f.orbit = 40; % many missing views!
%	f.orbit = 180; % use this for full angular sampling
	if 1 % parallel
		sg = sino_geom('par', 'nb', 140, 'na', 133, 'dr', ig.dx, ...
			'orbit', f.orbit, ... % missing views
			'orbit_start', (180-f.orbit)/2);
	else % fan
		sg = sino_geom('fan', 'nb', 250, 'na', 5, 'ds', ig.dx, ...
			'dso', 500, 'dod', 500, ...
			'dfs', inf, ... % flat fan
			'orbit', f.orbit, ... % missing views
			'orbit_start', (180-f.orbit)/2);
	end
	[sino t.pos] = ellipse_sino(sg, xell, 'oversample', 2);
	im(4, sg.s, sg.ad, sino, sprintf('sino: %d views', sg.na)), cbar
	if f.orbit < 180 && im
		axisy(0, 180)
		ytick([0 sg.orbit_start 90 sg.orbit+sg.orbit_start 180])
	end
prompt
end

if ~isvar('xfbp'), printm 'fbp'
	tmp = fbp2(sg, ig);
	xfbp = fbp2(sino, tmp);
	im(2, xfbp, 'fbp', clim), cbar
prompt
end

if ~isvar('G'), printm 'G'
	% informative support constraint:
	ig.mask = ellipse_im(ig, [0 0 95 125 0 1], 'oversample', 3) > 0;
	im(6, ig.mask, 'mask')

	if streq(sg.type, 'par')
		if has_mex_jf
			f.tab_type = {'square/strip', 'chat', 0, 'Ltab', ...
			1000, 'strip_width', ig.dx};
			G = Gtomo2_table(sg, ig, f.tab_type, 'nthread', 1);
		else
			G = Gtomo2_strip(sg, ig);
		end
	else
		G = Gtomo2_wtmex(sg, ig);
	end

	if 0
		t = G * xtrue;
		im(t), cbar
		max_percent_diff(t, sino)
	end
end

if ~isvar('R'), printm 'R'
	R = Reg1(ig.mask, 'type_denom', 'matlab', 'beta', 2^9, ...
		'pot_arg', {'hyper3', 1});
	psf = qpwls_psf(G, R, 1, ig.mask, 1);
%	xblur = conv2(xtrue, psf, 'same');
%	im(3, xblur, clim, 'xblur'), cbar
	im(3, psf, 'psf'), cbar
	clear xpcg xsps xiot
prompt
end

if ~isvar('xpcg'), printm 'xpcg'
	xinit = xfbp;
	f.niter = 10;

	P0 = 1;
	xpcg = pl_pcg_qs_ls(xinit(ig.mask), G, {sino(:), 1}, ...
		@wls_dercurv, R, 'precon', P0, ...
		'niter', f.niter, 'isave', f.niter);
	xpcg = ig.embed(xpcg);
	im(5, xpcg, clim, 'pwls-pcg'), cbar
prompt
end

if ~isvar('xsps'), printm 'xsps'
	if sg.na < 20
		Gb = G; % no subsets if very limited angle
	else
		Gb = Gblock(G, min(sg.na, 10));
	end
	xsps = pwls_sps_os(xinit(ig.mask), sino, ones(size(sino)), Gb, R, ...
		f.niter + 90*(sg.na < 20));
	xsps = ig.embed(xsps);
	im(6, xsps(:,:,end), clim, 'pwls-sps'), cbar
prompt
end

if ~isvar('xiot'), printm 'xiot'

	Gb = Gblock(G, min(sg.na, 26));
	xiot = pl_iot(xsps(ig.mask), Gb, {sino, ones(size(sino))}, R, ...
		'dercurv', 'wls', ...
		'riter', 2, ...
		'os', 5, ... % f.niter, ...
		'niter', f.niter, 'isave', 'last', ...
		'pixmin', 0, ...
		'chat', 0);
	xiot = ig.embed(xiot);

	im(4, xiot, clim, 'xiot'), cbar
prompt
end
