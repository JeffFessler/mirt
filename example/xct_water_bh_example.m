% xct_water_bh_example.m
%
% illustrate X-ray CT iterative reconstruction of an object
% consisting of water/air mixtures only (no bones) with simple
% beam-hardening corrections for polyenergetic spectrum.
%
% Copyright 2008-10-29, Jeff Fessler, University of Michigan

% todo: Poisson/OS version!

if ~isvar('ssino'), printm 'ssino'
	ig = image_geom('nx', 128, 'ny', 120, 'dx', 0.4);
	ig.mask = ig.circ > 0;
	xell = [0 0 15 10 0 1;
		-6 2 2 2 0 1.1;
		6 -2 2 2 0 0.9];
	[xtrue xell] = ellipse_im(ig, xell, 'oversample', 3);

	im plc 2 3, clim = [0 1.2];
	im(1, xtrue, clim, 'xtrue'), cbar

	sg = sino_geom('par', 'nb', 140, 'na', 100, 'dr', ig.dx);
	ssino = ellipse_sino(sg, xell, 'oversample', 2);
	im(4, sg.s, sg.ad, ssino, 's sino'), cbar
prompt
end


if ~isvar('ftab'), printm 'ftab'
	sls = de_ftab_sls('max', 50, 'n', 101);
	xrs = xray_read_spectra('poly1,140', ...
		'filters', {{'aluminum', 0.25, 'copper', 0.05}});
	mas = xray_read_mac('water');
	ftab = de_ftab(xrs, mas, 'sls', sls, 'ctype', 'newt', ...
		'ftype','exp', 'fit_args', {'kev', 10:5:160, 'mac', []});
%	im subplot 2
%	ftab.inv1.plot(ftab.fit);
prompt
end

if ~isvar('fsino'), printm 'fsino'
	fsino = ftab.fit.fmfun(ssino);
	im(5, sg.s, sg.ad, fsino, 'f sino'), cbar
prompt
end

if ~isvar('yi'), printm 'yi'
%	yi = fsino; % no noise
%	wi = 1; % unweighted

	f.I0 = 1e6; % high snr
	yi = poisson(f.I0 * exp(-fsino), 7) / f.I0;
	wi = yi;
	yi = -log(max(yi, 1/f.I0));
	im(6, sg.s, sg.ad, yi, 'yi sino'), cbar
	im(3, sg.s, sg.ad, wi, 'wi sino'), cbar
prompt
end

if ~isvar('fbpu'), printm 'fbpu'
	f.fbp = fbp2(sg, ig);
	tmp = yi / ftab.mac.bar; % scale correct for effective water mac
	fbpu = fbp2(tmp, f.fbp, 'window', 'hanning,0.8');
	im(2, fbpu, 'fbp uncorrected', clim), cbar
prompt
end

if ~isvar('fbpc'), printm 'fbpc'
	tmp = ftab.inv1.fun(yi); % water BH correction
	fbpc = fbp2(tmp, f.fbp, 'window', 'hanning,0.8');
	im(3, fbpc, 'fbp corrected', clim), cbar

	prompt
	if 1
		clf
		iy = ig.ny/2+1; ix = 1:ig.nx;
		pro = @(x) x(ix,iy);
		plot([pro(xtrue) pro(fbpc) pro(fbpu)])
		legend('true', 'FBP corrected', 'FBP uncorrected')
	prompt
	end
end

if ~isvar('G'), printm 'G'
	if has_mex_jf
		f.tab_type = {'square/strip', 'chat', 0, 'Ltab', ...
			1000, 'strip_width', ig.dx};
		G = Gtomo2_table(sg, ig, f.tab_type, 'nthread', 1);
	else
		G = Gtomo2_strip(sg, ig);
	end

	if 0
		tmp = G * xtrue;
		im(tmp), cbar
		max_percent_diff(tmp, ssino)
	end
end

if 0 % play with beta to find "desired" resolution
	f.beta = 2^1;
	R0 = Reg1(ig.mask, 'type_denom', 'matlab', 'beta', f.beta, ...
		'pot_arg', {'hyper3', 0.1});
	if 1
		psf = qpwls_psf(G, R0, 1, ig.mask, 1, 'loop', 1);
	end
return
end

if ~isvar('R'), printm 'R'
	% for simplicty use derivative at 0 for nonlinear term
%	tmp = wls_water_dercurv({0, 1, ftab, 1}, 0, [], 1, 1)
	tmp = ftab.fit.fgrad(0)
	ci = wi .* tmp^2; % adjust for effect of nonlinear BH function
	% Fessler 1996 spatial resolution properties...
	kappa = sqrt(div0(G' * ci, G' * sg.ones));
	R = Reg1(kappa, 'type_denom', 'matlab', 'beta', f.beta, ...
		'pot_arg', {'hyper3', 0.1});

	im(kappa)
prompt
end

%
% corrected iterative recon
%

if ~isvar('xc'), printm 'xc'
	f.niter = 400;
	xinit = fbpu;
%	xinit = fbpc > 0.9;
	data = {yi(:), wi(:), ftab, 1};
	dercurv = @wls_water_dercurv;
	tmp = pl_pcg_qs_ls(ig.maskit(xinit), G, ...
		data, dercurv, R, ...
		'niter', f.niter, 'isave', 'all');
	xc = ig.embed(tmp);

	clf, im(xc, clim)
	prompt
prompt
end

	if im
		clf, im pl 2 3
		im(1, xtrue, clim, 'true')
		im(2, fbpc, clim, 'fbp corrected')
		im(3, xc(:,:,end), clim, 'iterative')

		subplot(212)
		plot([pro(xtrue) pro(fbpc) pro(xc(:,:,end))])
		legend('true', 'FBP corrected', 'iterative')
	return
	end

%	movie2(xc, 'file', '/y/fessler/tmp3.avi')
