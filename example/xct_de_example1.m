% xct_de_example1.m
%
% illustrate dual energy (DE) X-ray CT reconstruction of an object
% consisting of water and bones, for two polyenergetic spectra.
% this is still a 'work in progress'
%
% Copyright 2008-12-19, Jeff Fessler, University of Michigan


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
	im plc 1 2, clim = [0 1.2];
	im(1, dens_true, clim, 'dens true'), cbar

	sg = sino_geom('par', 'nb', 140, 'na', 100, 'dr', ig.dx);
	wsino = ellipse_sino(sg, ellw, 'oversample', 2);
	bsino = ellipse_sino(sg, ellb, 'oversample', 2);

	ssino = cat(3, wsino, bsino);
	im(2, sg.s, sg.ad, ssino, 's sino'), cbar
prompt
end


if ~isvar('ftab'), printm 'ftab - generate x-ray spectra'
	xrs = xray_read_spectra('poly1,140,80', ...
		'filters', {{'aluminum', 0.25, 'copper', 0.05}});

	sls = de_ftab_sls('max', [50 40], 'n', [101 11]);
	mas = xray_read_mac({'water', 'bone'});
	ftab = de_ftab(xrs, mas, 'sls', sls, 'ctype', 'newt', ...
		'ctype', 'pre10', ...
		'ftype','exp', 'fit_args', {'kev', 10:5:160, 'mac', []});

	clf, ftab.plot_fm
prompt
end


if ~isvar('fsino'), printm 'fsino'
	fsino = ftab.fit.fmfun(ssino);
	im(sg.s, sg.ad, fsino, 'f sino'), cbar
prompt
end


if ~isvar('yi'), printm 'yi'
%	yi = fsino; % no noise
%	wi = 1; % unweighted

	f.I0 = 1e6; % high snr for now
	f.I0 = f.I0 * xrs.I / xrs.I(1); % I0 proportional to spectrum integral
	yi = zeros(size(ssino));
	wi = zeros(size(ssino));
	for mm=1:size(ssino, 3)
		tmp = poisson(f.I0(mm) * exp(-fsino(:,:,mm)), 7) / f.I0(mm);
		wi(:,:,mm) = tmp;
		yi(:,:,mm) = -log(max(tmp, 1/f.I0(mm)));
	end
	im plc 1 2
	im(1, sg.s, sg.ad, yi, 'yi sino'), cbar
	im(2, sg.s, sg.ad, wi, 'wi sino'), cbar
prompt
end


if 0 % picture of fast kvp switching
	im pl 1 2
	tmp = fsino(1:1:end,1:2:end,:);
	tmp(:,2:2:end,1) = 0;
	tmp(:,1:2:end,2) = 0;
	im(1, tmp(:,:,1), 'Low Energy Sinogram')
	axis normal, xtick off, ytick off
	im(2, tmp(:,:,2), 'High Energy Sinogram')
	axis normal, xtick off, ytick off
return
end


if ~isvar('fbpc'), printm 'fbpc' % water corrected only
	tmp = ftab.inv1.fun(yi); % water BH correction
	f.fbp = fbp2(sg, ig);
	fbpc = fbp2(tmp, f.fbp, 'window', 'hanning,0.8');
	im clf, im(fbpc, 'fbp corrected', clim), cbar

	prompt
	if 1
		clf
		iy = ig.ny/2+1; ix = 1:ig.nx;
		pro = @(x) x(ix,iy); % profile
		% pseudo-density
		xtrue = wtrue + btrue * ftab.mac.bar(2) / ftab.mac.bar(1);
		plot([pro(xtrue) pro(dens_true) pro(fbpc)])
		legend('xtrue', 'dens true', 'FBP corrected', 'FBP uncorrected')
	prompt
	end
end


if ~isvar('shat'), printm 's hat : sinogram domain DE decomposition'
	shat = ftab.inv2.fun(yi);
	im(shat), cbar
prompt
end


if ~isvar('fbp'), printm 'fbp DE'
	fbp = fbp2(shat, f.fbp, 'window', 'hanning,0.8');
	im clf, im(fbp, 'fbp DE', clim), cbar
prompt
end


% iterative recon for DE CT covered by US Patent 6,754,298
if exist('xct_de_example1_p2.m') == 2 % mfile for UM only
	xct_de_example1_p2
end
