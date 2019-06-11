% ir_denoise_admm1_test
% illustrate using "ADMM" for image denoising (aka Split Bregman algorithm)

if ~isvar('yi'), printm 'generate noisy image'
%	xt = single(imread('cameraman.tif')');
%	sig = 6; clim = [0 250];
	xt = ellipse_im(2^7);
	sig = 2^-1; clim = [0 9];

	snr = @(x) -20*log(norm(x(:)-xt(:)) / norm(xt));

	rng(7)
	yi = xt + sig * randn(size(xt));
	im plc 2 2
	im(1, xt, clim, 'True')
	im(2, yi, clim, 'Noisy')
	xlabelf('SNR = %4.1f dB', snr(yi))
end


if ~isvar('xs'), printm 'run admm'
	[xs snrs] = ir_denoise_admm1(yi, 'beta', 2^-1, ...
		'stop_diff_tol', 1e-4, 'chat', 0, ...
		'isave', 'all', ...
		'userfun', @(x, iter) snr(x), ...
		'niter', 2^5, 'shrink', [], 'rho', 1);
	xh = xs(:,:,end);
end


if 1 % figures
%	im plc 1 3
%	im(1, xt, clim)
%	im(2, yi, clim)
	im(3, xh, clim, 'Denoised')
	xlabelf('SNR = %4.1f dB', snr(xh))

	im subplot 4
	iters = 0:numel(snrs)-1;
	plot(iters, snrs, '-o')
	axis([0 20 20 55])
	xlabel 'ADMM iteration'
	ylabel 'SNR (dB)'
	ytick([25 53])
%	ir_savefig cw ir_denoise_admm1_test
end


if 0 % movie of iterates
	tmp = xs(:,:,[1:20 end]);
	niter = size(tmp,3);
	bar = size(tmp,2) / niter;
	for ii=1:niter
                tmp(1:3, floor([1:bar]+(ii-1)*bar), ii) = 8;
        end
	clf, im(tmp)
%	movie2(tmp, 'clim', clim, 'fps', 10)
end
