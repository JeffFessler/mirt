% xct_bh_example.m
%
% illustrate X-ray CT iterative reconstruction of an object
% consisting of water and bones, for a polyenergetic spectrum.
%
% Copyright 2008-12-05, Jeff Fessler, University of Michigan

warn 'work in progress'

if ~isvar('ssino'), printm 'ssino'
	ig = image_geom('nx', 128, 'ny', 120, 'dx', 0.4);
	ig.mask = ig.circ > 0;
	ellw = [0 0 15 10 0 1;
		-6 0 2 2 0 -1;
		6 -0 2 2 0 -1];
	wtrue = ellipse_im(ig, ellw, 'oversample', 3);

	ellb = [-6 0 2 2 0 1.9;
		6 -0 2 2 0 1.9];
	btrue = ellipse_im(ig, ellb, 'oversample', 3);

	dens_true = wtrue + btrue;
	im plc 2 3, clim = [0 1.2];
	im(1, dens_true, clim, 'dens true'), cbar

	sg = sino_geom('par', 'nb', 140, 'na', 100, 'dr', ig.dx);
	wsino = ellipse_sino(sg, ellw, 'oversample', 2);
	bsino = ellipse_sino(sg, ellb, 'oversample', 2);

	ssino = cat(3, wsino, bsino);
	im(4, sg.s, sg.ad, ssino, 's sino'), cbar
prompt
end


if ~isvar('ftab2'), printm 'ftab2'
	xrs = xray_read_spectra('poly1,140', ...
		'filters', {{'aluminum', 0.25, 'copper', 0.05}});

	sls1 = de_ftab_sls('max', 50, 'n', 101);
	mas1 = xray_read_mac('water');
	ftab1 = de_ftab(xrs, mas1, 'sls', sls1, 'ctype', 'newt', ...
		'ftype','exp', 'fit_args', {'kev', 10:5:160, 'mac', []});

	sls2 = de_ftab_sls('max', [50 40], 'n', [101 11]);
	mas2 = xray_read_mac({'water', 'bone'});
	ftab2 = de_ftab(xrs, mas2, 'sls', sls2, 'ctype', 'newt', ...
		'ftype','exp', 'fit_args', {'kev', 10:5:160, 'mac', []});

	ftab2.plot_fm
%	ftab1.inv1.plot % todo: fails - fix!
	ftab1.inv1.plot(ftab1.fit);
%	im subplot 2
prompt
end

if 0
	clf
	tmp = ftab1.inv1.fun(ftab2.fm);
	plot(sls1.sl{1}, tmp)
	xlabel 's1', ylabel ''
return
end

if ~isvar('fsino'), printm 'fsino'
	fsino = ftab2.fit.fmfun(ssino);
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
	tmp = yi / ftab2.mac.bar(1); % scale correct for effective water mac
	fbpu = fbp2(tmp, f.fbp, 'window', 'hanning,0.8');
	im(2, fbpu, 'fbp uncorrected', clim), cbar
prompt
end

if ~isvar('fbpc'), printm 'fbpc'
	tmp = ftab1.inv1.fun(yi); % water BH correction
	fbpc = fbp2(tmp, f.fbp, 'window', 'hanning,0.8');
	im(3, fbpc, 'fbp corrected', clim), cbar

	prompt
	if 1
		clf
		iy = ig.ny/2+1; ix = 1:ig.nx;
		pro = @(x) x(ix,iy);
		% pseudo-density
		xtrue = wtrue + btrue * ftab2.mac.bar(2) / ftab2.mac.bar(1);
		plot([pro(xtrue) pro(dens_true) pro(fbpc) pro(fbpu)])
		legend('xtrue', 'dens true', 'FBP corrected', 'FBP uncorrected')
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

if ~isvar('fbpb'), printm 'fbpb' % "bone only recon"
	tmp = fbpc .* (fbpc > 1.3); % bone only image
	im(tmp)

	tmp = G * tmp;
	tmp = tmp.^2;
	fbpb = fbp2(tmp, f.fbp, 'window', 'hanning,0.8');
	im(fbpb)
return
end

if 1 % test it - it kind of works!
	scale =  0.005; % empirical value
	fbp1 = fbpc + scale * fbpb;

	if 1
		clf
		iy = ig.ny/2+1; ix = 1:ig.nx;
		pro = @(x) x(ix,iy);
		% pseudo-density
		xtrue = wtrue + btrue * ftab2.mac.bar(2) / ftab2.mac.bar(1);
		plot([pro(xtrue) pro(dens_true) pro(fbpc) pro(fbp1)])
		legend('xtrue', 'dens true', 'FBP corrected', 'FBP1 corrected')
	end

	clim = 1 + [-1 1] * 0.10;
	im(stackup(xtrue, fbpu, fbpc, fbp1), clim), cbar
return
end


if ~isvar('R'), printm 'R'
	% Fessler 1996 spatial resolution properties...
	kappa = sqrt(div0(G' * wi, G' * sg.ones));
	R = Reg1(kappa, 'type_denom', 'matlab', 'beta', 2^3, ...
		'pot_arg', {'hyper3', 0.1});
	if 1
		qpwls_psf(G, R, 1, ig.mask, Gdiag(wi));
	end
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
