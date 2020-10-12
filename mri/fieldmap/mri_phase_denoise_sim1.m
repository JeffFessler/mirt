% mri_phase_denoise_sim1.m
% simulation using mri_phase_denoise()

if ~isvar('xtrue'), printm 'xtrue'
 if 0
	% first form "true" smooth phase map
	f.pmap = './brain,pmap,o2,l2b-6.fld';
	if exist(f.pmap, 'file')
		xtrue = fld_read(f.pmap);
	else
		xtrue = mri_phase_denoise(yi, 'pl', 1, 'order', 2, 'l2b', -6);
		im(xtrue, 'PL-CG phase, o=2'), cbar
%		fld_write(f.pmap, xtrue);
	end
	ig = image_geom('nx', 64, 'ny', 64, 'dx', 1);
	ymask = ellipse_im(ig, [0 0 27 35 0 1], 'oversample', 3) > 0;
	yi_lim = [0 max(abs(yi(:)))];

 else % generate simulated field map and image
%	xtrue = -20 * ones(64,64);
	xtrue = zeros(64,64);
%	xtrue(33+[-15:15], 20:35) = 90;
	xtrue(33+[-13:13], 20:38) = 1;
	h = ones(5); h = h/sum(h(:));
	xtrue = conv2(xtrue, h, 'same');
%	xtrue = -20 + 110 * xtrue;
	xtrue = -20 + 50 * xtrue;
	xtrue(28, 47) = max(xtrue(:));
	f.etd = 2e-3; % echo-time difference
	xtrue = xtrue * 2*pi * f.etd; % phase in radians
	f.dir = [path_find_dir('mri') '/../data/mri/'];
	f.xtrue = [f.dir 'brainweb_t1.jpg'];
	yi = double6(imread([f.dir 'brainweb_t1.jpg']))';
	yi = yi(2:end-1, 2:end-1);
	yi = downsample2(yi, 4);
	ig = image_geom('nx', 64, 'ny', 64, 'dx', 1);
	ymask = ellipse_im(ig, [0 0 21 25 0 1], 'oversample', 3) > 0;
	yi_lim = [0 180];
 end
end

if ~isvar('xo'), printm 'xo'
	rng(0)
	ytrue = abs(yi) .* exp(1i * xtrue);
	yi = ytrue + 5 * (randn(size(ytrue)) + 1i * randn(size(ytrue)));
	printm('snr=%g dB', 10*log10(rms(yi(:)) / rms(yi(:)-ytrue(:))))

	im plc 3 4; cbar off, im tickoff
	scale = 1 / (2*pi*f.etd); % 2ms echo time difference
%	clim = [-40 120];
	clim = round(minmax(xtrue)'/f.etd/2/pi + [-5 5]);
	elim = [-20 20];
	im(1, scale*xtrue .* ymask, 'True field map', clim), cbar
	t = (conv2(single(ymask), ones(3), 'same') > 0) - ymask;
	im(5, abs(yi) + 0*800*t, '|yi|', yi_lim), cbar

	xo = angle(yi);
	xo = xo .* ymask;
	nrms(xo(ymask), xtrue(ymask))
	err = scale * (xo-xtrue);
	printm('conv rms = %g, max = %g', rms(err(ymask)), max(abs(err(ymask))))

	im(2, scale*xo, 'Conventional', clim), cbar
	im(6, err.*ymask, elim, 'error'), cbar
	xlabelf('RMSE = %4.1f Hz', rms(err(ymask)))

	if 0 % for movie
		xall = mri_phase_denoise(yi, 'pl', 1, 'order', 2, 'l2b', -6, 'isave', 0:150);
		fld_write('xall', xall)
		keyboard
	end
end

if ~isvar('xpwls'), % run PL and PWLS estimator
	cpu etic
	xpl = mri_phase_denoise(yi, 'pl', 1, 'order', 2, 'l2b', -6);
	cpu etoc 'PL time'
	xpl = xpl .* ymask;
	nrms(xpl(ymask), xtrue(ymask))
	err = scale * (xpl-xtrue);
	printm('pl rms = %g, max = %g', rms(err(ymask)), max(abs(err(ymask))))

	im(3, scale*xpl, 'PL', clim), cbar
	im(7, err.*ymask, elim, 'error'), cbar
	xlabelf('RMSE = %4.1f Hz', rms(err(ymask)))

	cpu etic
	xpwls = mri_phase_denoise(yi, 'wi_ml', 1, 'order', 2, 'l2b', -6);
	cpu etoc 'PWLS time'
	xpwls = xpwls .* ymask;
	nrms(xpwls(ymask), xtrue(ymask))
	err = scale * (xpwls-xtrue);
	printm('pwls rms = %g, max = %g', rms(err(ymask)), max(abs(err(ymask))))

	im(4, scale*xpwls, 'PWLS (ML)', clim), cbar
	im(8, scale*(xpwls-xtrue.*ymask), elim, 'error'), cbar
	xlabelf('RMSE = %4.1f Hz', rms(err(ymask)))
prompt
end

return

if 0
	xptls = mri_phase_denoise(yi, 'wi_ml', 0, 'order', 2, 'l2b', -6);
	xptls = xptls .* ymask;
	nrms(xptls(ymask), xtrue(ymask))
	err = scale * (xptls-xtrue);
	printm('ptls rms = %g, max = %g', rms(err(ymask)), max(abs(err(ymask))))

	im(5, scale*xptls, 'PWLS (simple)', clim), cbar
	subplot(3,5,10)
	im(scale*(xptls-xtrue.*ymask), elim, 'error'), cbar
end

%	ir_savefig fig_mri_phase_denoise_sim1
%	ir_savefig fig_mri_phase_denoise_sim2
%keyboard

% examine effect on EPI recon

if ~isvar('f0fft'), printm 'f0fft'
	ftrue = zeros(size(xtrue));
	ftrue(7:1:end-7,7:1:end-7) = 1;
	ftrue(7:5:end-7,7:1:end-7) = 2;
	ftrue(7:1:end-7,7:5:end-7) = 2;
	[nx ny] = size(ftrue);
	[kx ky] = ndgrid([-nx/2:nx/2-1],[-ny/2:ny/2-1]);
	mask = true(nx,ny);
	ti = [0:(nx*ny-1)]/nx/ny * 30e-3; % readout
	arg = {[kx(:) ky(:)], mask, 'ti', ti, 'L', 8, 'zmap'};

	Gtrue = Gmri(arg{:}, 1i * xtrue / f.etd);

	st = Gtrue * ftrue(:); % signal
	st = reshape(st, size(ftrue)); % EPI

	f0fft = fftshift(ifft2(fftshift(st))); % uncorrected recon
	clf, im([ftrue; f0fft]), cbar

	clear Go
prompt
end

if ~isvar('Go'), printm 'Go'
	G0 = Gmri(arg{:}, 0i * xtrue / f.etd);
%	T0 = build_gram(G0, 1);
	Gpl = Gmri(arg{:}, 1i * xpl / f.etd);
%	Tpl = build_gram(Gpl, 1);
	Gpwls = Gmri(arg{:}, 1i * xpwls / f.etd);
%	Tpwls = build_gram(Gpwls, 1);
	Go = Gmri(arg{:}, 1i * xo / f.etd);
%	To = build_gram(Go, 1);

	clear fo
end

if ~isvar('R') || 0
	R = Robject(mask, 'beta', 2^7, 'order', 2);
	qpwls_psf(G0, R, 1, mask);
	clear fo
prompt
end

if ~isvar('fo'), printm 'fo1'
	f.niter = 20;

	init = f0fft(mask);
	fpl = qpwls_pcg1(init, Gpl, 1, st(:), R.C, 'niter', f.niter);
	fpl = embed(fpl, mask);
	im(fpl), cbar
	fpwls = qpwls_pcg1(init, Gpwls, 1, st(:), R.C, 'niter', f.niter);
	fpwls = embed(fpwls, mask);
	im(fpwls), cbar
	fo = qpwls_pcg1(init, Go, 1, st(:), R.C, 'niter', f.niter);
	fo = embed(fo, mask);
	im(fo), cbar
	f0 = qpwls_pcg1(init, G0, 1, st(:), R.C, 'niter', f.niter);
	f0 = embed(f0, mask);
	ft = qpwls_pcg1(init, Gtrue, 1, st(:), R.C, 'niter', f.niter);
	ft = embed(ft, mask);
prompt
end

if ~isvar('fo'), printm 'fo2'
	f.niter = 20;

	init = f0fft(mask);
	fpl = qpwls_pcg2(init, Tpl, Gpl'*st(:), R.C, 'niter', f.niter);
	fpl = embed(fpl, mask);
	im(fpl), cbar
	fpwls = qpwls_pcg2(init, To, Gpwls'*st(:), R.C, 'niter', f.niter);
	fpwls = embed(fpwls, mask);
	im(fpwls), cbar
	fo = qpwls_pcg2(init, To, Go'*st(:), R.C, 'niter', f.niter);
	fo = embed(fo, mask);
	im(fo), cbar
	f0 = qpwls_pcg2(init, T0, G0'*st(:), R.C, 'niter', f.niter);
	f0 = embed(f0, mask);
end

if 1
	im plc 3 5
	flim = [0 2.2];
%	im(1, 0*xtrue, clim, 'None')
	im(1, ftrue, ' ')
	im(2, xtrue*scale, clim, 'True fieldmap'), cbar below
	im(3, xo*scale, clim, 'Conventional')
	im(4, xpl*scale, clim, 'PL')
	im(5, xpwls*scale, clim, 'PWLS')

	xl = @(ft,f0) xlabelf('NRMS %4.2f%%', nrms(ft, f0));
%	xl = @(x1, x2) sum(col(x1+x2)); % dummy to hide
	im(6, abs(f0), flim, 'None')
	xl(ftrue, f0);
	im(7, abs(ft), flim, ' ')
	xl(ftrue, ft);
	%im(6, f0fft, 'f0fft', flim)
	%im(5, ftrue, flim)
	im(8, abs(fo), flim, ' ')
	xl(ftrue, fo);
	im(9, abs(fpl), flim, ' ')
	xl(ftrue, fpl);
	im(10, abs(fpwls), flim, ' ')
	xl(ftrue, fpwls);
%	ir_savefig mri_phase_denoise_sim2_epi
end
