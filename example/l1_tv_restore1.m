% l1_tv_restore1
% illustrate image restoration using a l1-type data fit term
% (for noise robustness) and a (anisotropic) TV-type regularizer
% based on l1_regress_example.m
% Copyright 2010-05-25, Jeff Fessler, University of Michigan

if ~isvar('xtrue')
	f.dir = [path_find_dir('mri') '/../data/mri/'];
	f.xtrue = [f.dir 'brainweb_t1.jpg'];
	xtrue = single(imread(f.xtrue)');
	xtrue = xtrue(2:end-1,2:end-1);
	im(xtrue), cbar
end

	snr = @(a,b) 20 * log10(norm(a(:))/norm(b(:)));

if ~isvar('A')
%	psf = fspecial('Gaussian', [7 7], 5)
	psf = gaussian_kernel(5 * sqrt(log(256)), 3);
	psf = psf * psf'; psf = psf / sum(psf(:));
	mask = true(size(xtrue));
	A = Gblur(mask, 'psf', psf);

	yb = A * xtrue;
	im(yb), cbar
	f.snr_db = 25;
	sig = 10^(-f.snr_db/20) * norm(yb(:)) / sqrt(numel(yb));

	rng(0)
	yi = yb + sig * randn(size(yb));
	printm('data snr = %g dB', snr(yb, yi-yb))

	im(yi), cbar
end


if ~isvar('xh')
	xh = l1_tv_restore1_fun(yi, A, 'l2b', -2);
	xh = embed(xh, mask);
	im(xh), cbar
end

	printm('raw snr = %g dB', snr(xtrue, yi - xtrue))
	printm('new snr = %g dB', snr(xtrue, xh - xtrue))

if 0
	beta = 1e-3;
	gamma = 1e-6;
	mu = 1/25; % for 25% random-valued noise
%	sqrt((y - A f)^2 + gamma) + mu * sqrt((Cx f)^2 + (Cy *f)^2 + beta)
end
