% ir_denoise_test_sep1
% Examine denoising with ordinary vs separable finite differences
% i.e. |Ch X|_1 + |Cv X|_1 vs |Ch X Cv|_1
% Uses ADMM algorithm to optimize denoising cost function.

if ~isvar('yi'), printm 'generate noisy image'
	xt = ellipse_im(2^7);
	xt = zeros(2^6); xt(20:60,20:48) = 9;
	sig = 2^0; clim = [0 9];

	snr = @(x) -20*log(norm(x(:)-xt(:)) / norm(xt(:)));

	rng(7)
	yi = xt + sig * randn(size(xt));

	im plc 2 4
	im(1, xt, clim, 'True')
	im(2, yi, clim, 'Noisy')
	xlabelf('SNR = %4.1f dB', snr(yi))
end


if ~isvar('Cs'), printm 'C2 and Cs'
	arg.Ch = Cdiffs(size(xt,1), 'type_diff', 'circshift', 'offsets', 1);
	arg.Cv = Cdiffs(size(xt,2), 'type_diff', 'circshift', 'offsets', 1);
	Cs = fatrix2('arg', arg, 'idim', size(xt), 'odim', size(xt), ...
		'forw', @(arg, x) arg.Ch * x * arg.Cv, ... % separable!
		'back', @(arg, y) arg.Ch' * y * arg.Cv'); 
	C2 = Cdiffs(size(xt), 'type_diff', 'circshift', 'offsets', '2d:hv');

	im(3, 'col', 1, C2 * xt, 'Usual'), cbar
	im(4, Cs * xt, 'Separable'), cbar
	drawnow
end


if ~isvar('x2'), printm 'run admm with usual 2d finite differences'
	[x2s snr2] = ir_denoise_admm1(yi, 'beta', 2^+1, ...
		'stop_diff_tol', 1e-4, 'chat', 0, ...
		'isave', 'all', ...
		'userfun', @(x, iter) snr(x), ...
		'niter', 2^5, 'shrink', [], 'rho', 1);
	x2 = x2s(:,:,end);
end

if ~isvar('x1'), printm 'run admm with separable finite differences'
	[x1s snr1] = ir_denoise_admm1(yi, 'beta', 2^+1, ...
		'stop_diff_tol', 1e-4, 'chat', 0, ...
		'C', 'sep1', ...
		'isave', 'all', ...
		'userfun', @(x, iter) snr(x), ...
		'niter', 2^5, 'shrink', [], 'rho', 1);
	x1 = x1s(:,:,end);
end


if 1 % figures
	im(7, x2, clim, 'Denoised usual')
	xlabelf('SNR = %4.1f dB', snr(x2))
	im(8, x1, clim, 'Denoised separ.')
	xlabelf('SNR = %4.1f dB', snr(x1))

%	im subplot 5
	subplot(223)
	iters = 0:numel(snr1)-1;
	plot(iters, snr2, '-o', iters, snr1, '-x')
	snr_lim = [floor(min([snr1; snr2])) ceil(max([snr1; snr2]))];
	axis([0 numel(iters) snr_lim])
	ytick(snr_lim)
	xlabelf 'ADMM iteration'
	ylabelf 'SNR (dB)'
	legend('usual', 'separable', 'location', 'southeast')
end
