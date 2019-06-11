% restore_example2.m
% Example of edge-preserving image restoration
% using Gblur object for system model (shift-invariant blur)
% and Reg1 for nonquadratic regularization.
%
% Copyright 2005-5-18, Jeff Fessler, University of Michigan

if ~isvar('xtrue'), printm 'read image'
	ig = image_geom('nx', 2^8, 'dx', 1);
	xtrue = 256/9*ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
end

% PSF from figueiredo:05:abo
if ~isvar('psf'), printm 'psf'
	psf = -7:7;
	psf = ndgrid(-7:7, -7:7);
	psf = 1 ./ (1 + psf.^2 + (psf').^2);
	psf = psf / sum(psf(:)); % normalize to unity DC response
	im plc 2 2
	im(1, psf, 'psf')
end

if ~isvar('A'), printm 'system model'
	A = Gblur(ig.mask, 'psf', psf);
end

if ~isvar('yi')
	y0 = conv2(xtrue, psf, 'same');

	rng(0)
	estd = sqrt(2);
	yi = y0 + estd * randn(size(y0));

	clim = [0 255];
	im(2, yi, 'yi', clim)
end

if ~isvar('xnpls') || 1
	f.l2b_n = -2; % hyper3
	f.delta = 0.3;
	Rn = Reg1(ig.mask, 'type_denom', 'matlab', ...
		'pot_arg', {'hyper3', f.delta}, 'beta', 2^f.l2b_n);

	f.niter = 60;
	xinit = yi(ig.mask);
	xnpls = pwls_sps_os(xinit, yi(:), [], A, Rn, ...
		f.niter, [0 inf], [], [], 1);

	snr_improve = 10 * log10(sum(col(yi-xtrue).^2) ./ ...
		sum((xnpls-xtrue(:)*ones(1,f.niter)).^2));
	if im
		subplot(224), plot(0:f.niter-1, snr_improve, '-o')
		xlabel 'iter', ylabel 'SNR improvement [dB]'
	end

	xnpls = embed(xnpls, ig.mask);
	im(3, xnpls(:,:,end), 'xnpls', clim)
end
