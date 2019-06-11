% restore_example.m
% Example of edge-preserving image restoration
% using Gblur object for system model (shift-invariant blur)
% and nonquadratic regularization using Reg1.
%
% Copyright 2002-5-10, Jeff Fessler, University of Michigan

if ~isvar('xtrue'), printm 'read image'
%	xtrue = double(imread('fig-tomos-old.gif'))';
%	xtrue = (1 - xtrue / 3) * 100;
%	xtrue = xtrue(50+[1:1100], 60+[1:360], :);
	try
		xtrue = double(imread('fig,tomos.gif'))';
	catch
		xtrue = double(imread('fig,tomos.tif'))';
	end
	xtrue = (1 - xtrue / max(xtrue(:))) * 100;
	xtrue = xtrue(50+[1:1000], 0+[1:400], :);
	kern = gaussian_kernel(4);
	xtrue = conv2(xtrue, kern, 'same');
	xtrue = conv2(xtrue, kern', 'same');
	xtrue = xtrue(2:4:end,2:4:end);	% downsample for speed
%	xtrue = xtrue(1:2:end,:) + xtrue(2:2:end,:);
%	xtrue = xtrue(:,1:2:end) + xtrue (:,2:2:end);
	[nx,ny] = size(xtrue);
end


% this PSF makes visible blur yet is still invertible (for illustration)
if ~isvar('psf'), printm 'psf'
	psf = ones(5,1);
	psf((end+1)/2) = 5;
	psf = psf / sum(psf(:));
	psf = psf * psf';

	psfpad = zeros(nx,ny);
	psfpad(1:size(psf,1), 1:size(psf,2)) = psf;
	psfpad = psfpad([3:nx 1 2], :);
	psfpad = psfpad(:, [3:ny 1 2]);
	Psf = reale(fft2(psfpad));
end

clim = [0 100];

if ~isvar('xinv')
	y0 = conv2(xtrue, psf, 'same');

	x0 = reale(ifft2((1 ./ Psf) .* fft2(y0)));
	x0 = max(x0,0);

	rng(0)
	estd = 10;
	yi = y0 + estd * randn(size(y0));

	xinv = reale(ifft2((1 ./ Psf) .* fft2(yi)));

	xml = 100 * (xinv > 50);
end

if ~isvar('G')
	G = Gblur(true(nx, ny), 'psf', psf);
	mask = G.arg.mask;
end

if ~isvar('xqpls')
	f.l2b_q = -1;
	f.niter = 20;
	Rq = Reg1(mask, 'type_denom', 'matlab', ...
		'pot_arg', {'quad'}, 'beta', 2^f.l2b_q);

	xinit = yi;

	if 0 % timing new vs old
		Ro = Robject(mask, 'type_denom', 'matlab', ...
			'potential', 'quad', 'beta', 2^f.l2b_q);
		% warm up
		pwls_sps_os(xinit(:), yi(:), [], G, Rq, 1);
		pwls_sps_os(xinit(:), yi(:), [], G, Ro, 1);
		cpu etic
		xo = pwls_sps_os(xinit(:), yi(:), [], G, Ro, f.niter);
		cpu etoc old
		xo = embed(xo, mask);

		cpu etic
		xqpls = pwls_sps_os(xinit(:), yi(:), [], G, Rq, f.niter);
		cpu etoc new
		xqpls = embed(xqpls, mask);
		equivs(xqpls, xo) % almost same except at edges
	return
	end

	xqpls = pwls_sps_os(xinit(:), yi(:), [], G, Rq, f.niter);
	xqpls = embed(xqpls(:,end), mask);
end

if ~isvar('xnpls')
%	f.l2b_n = 1.5; % cauchy
	f.l2b_n = 8;
	Rn = Reg1(mask, 'type_denom', 'matlab', ...
		'pot_arg', {'hyper3', 0.1}, 'beta', 2^f.l2b_n);
%		'potential', 'cauchy', 'delta', 10);
	xinit = yi;
%	cpu etic
	xnpls = pwls_sps_os(xinit(:), yi(:), [], G, Rn, 3*f.niter);
%	cpu etoc new
	xnpls = embed(xnpls(:,end), mask);
end

if ~isvar('xo') && 0 % test old vs new
	Ro = Robject(mask, 'type_denom', 'matlab', ...
		'potential', 'hyper3', 'beta', 2^f.l2b_n, 'delta', 0.1);
	cpu etic
	xo = pwls_sps_os(xinit(:), yi(:), [], G, Ro, 3*f.niter);
	cpu etoc old
	xo = embed(xo(:,end), mask);
	max_percent_diff(xnpls, xo)
return
end

% wiener filter
if ~isvar('xwien'), printm 'wiener filter'
	yi0 = yi;
%	yi0 = eye(4);
	yi0 = yi0 - mean(yi0(:)); % zero-mean
	acorr2 = @(x) conv2(x, rot90(conj(x), 2)); % from xcorr2
	Ry = acorr2(yi0) / numel(yi0); % empirical autocorrelation
	[mx my] = size(yi0);
	Ry = Ry([-mx/2:mx/2-1]+mx,[-my/2:my/2-1]+my);
	Ry(:,1) = 0; Ry(1,:) = 0;
	Py = reale(fftshift(fft2(ifftshift(Ry))), 'warn');
	Px = max(Py-estd^2, 0);
	H = Px ./ (Px + estd^2); % Wiener filter
	xwien = reale(ifft2(fft2(yi) .* fftshift(H)));
%	imax2(Ry);
%	plot(Ry(:,(end+1)/2))
%	im(max(Ry-estd^2, 0))
prompt
end

% individual figures for preface example
if 0
	clf
	set(0, 'DefaultAxesFontSize', 22)
	im('notick', xtrue, clim), colormap(1-gray)
	title 'True image x'
	ir_savefig 'fig_front_xtrue'

	im('notick', yi, clim), colormap(1-gray)
	title 'Noisy/blurry data y'
	ir_savefig 'fig_front_y'

	im('notick', xinv, clim), colormap(1-gray)
	title 'Inverse filter estimated x'
	ir_savefig 'fig_front_xinv'

	im('notick', xnpls, clim), colormap(1-gray)
	title 'Regularized, statistical estimated x'
	ir_savefig 'fig_front_xhat'
return
end

ylab = @ (xh) ylabelf('%5.2f', nrms(xh(:), xtrue(:)));

if 1 && im
	im plc 5 2
	set(0, 'DefaultAxesFontSize', 16)
	im(1, 'notick', xtrue, clim), ylab(xtrue)
	title 'True x'

	im(2, 'notick', y0, clim), ylab(y0)
	title 'Blurry y (no noise)'

	im(3, 'notick', x0, clim), ylab(x0)
	title 'Inverse, no noise'

	im(4, 'notick', yi, clim), ylab(yi)
	title 'Noisy/blurry y'

	im(5, 'notick', xinv, clim), ylab(xinv)
	title 'Inverse, with noise'

	im(6, 'notick', xml, clim), ylab(xml)
	title 'Thresholded Inverse'

	im(7, 'notick', xqpls, clim), ylab(xqpls)
	title 'Quadratic Regularized LS'

	im(8, 'notick', xnpls, clim), ylab(xnpls)
	title 'Non-Quad. Reg. LS'

	im(9, 'notick', xwien, clim), ylab(xwien)
	title 'Wiener Filter'

	colormap(1-gray)

%	set(gcf, 'papersize', [8.5 8]), orient tall
%	ir_savefig 'fig_tomos'
return
end

if 1 && im
	c = (Rn.C1) * xnpls(mask);
	c = reshape(c, [size(mask) 4]);
	c = Rn.wpot(1, c);
	im plc 3 2
	im(1, 'notick', c(:,:,1), 'Horizontal')
	im(2, 'notick', c(:,:,2), 'Vertical')
	im(3, 'notick', c(:,:,3), 'Diagonal Up-Left')
	im(4, 'notick', c(:,:,4), 'Diagonal Up-Right')
%	set(gcf, 'papersize', [8.5 6]), orient tall
%	ir_savefig 'fig_tomos_wpot'
end
