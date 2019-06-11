%| Example regularized estimation of MR coil sensitivity maps.
%| mri_sensemap_demo1.m
%|
%| This is an older slower method that is made obsolete by the new version:
%| ir_mri_sensemap_admm.m
%| (However, for small 2D problems the cholesky approach below is best.)
%|
%| Copyright 2005-9-12, Jeff Fessler, University of Michigan
%| 2018-03-01 JF added 'ichol' test

%% true object image
if ~isvar('ftrue')
	f.dir = [path_find_dir('mri') '/../data/mri/'];
	f.xtrue = [f.dir 'brainweb_t1.jpg'];
	ftrue = double(imread(f.xtrue)');
	ftrue = ftrue(2:end-1,2:end-1); % make it 256^2
	ftrue = downsample2(ftrue, 4); % now 64^2
	[nx, ny] = size(ftrue);
	atrue = 2*pi * (-0.5+([1:nx]'/nx * [1:ny]/ny).^0.5); % smooth phase
%	im('hsv', atrue, [-pi pi]), colorbar
	ftrue = ftrue .* exp(1i * atrue); % phase
%	im('hsv', angle(ftrue))

	if 0 % simple synthetic rho map
		ftrue = zeros(size(ftrue));
		ftrue(nx/4:3*nx/4,ny/4:3*ny/4) = 2^7;
		ftrue(5+nx/2+[-2:2],ny/2+[-10:10]) = 2^2;
	%	ftrue = 128 * ones(size(ftrue));
	end

	im plc 3 3
	im(1, abs(ftrue), '|f| true'), cbar
	drawnow
end


%% noisy body and surface coil images
if ~isvar('zlj'), printm 'zlj'
	im row 2
	ncoil = 4;
	smap = ir_mri_sensemap_sim('nx', nx, 'ny', ny, 'dx', 192/nx, ...
		'ncoil', ncoil, 'orbit', 90*ncoil, 'rcoil', 100);
	if ncoil == 4, smap = smap(:,:,[4 1 3 2]); end
%	smap(nx/2+[-1:1]*10,ny/2+[-1:1]*10,1) = 1.5;

	f.sig = 1.0; % noise in image: bodycoil and array coils
	rng(0)
	yj = ftrue + f.sig * (randn(size(ftrue)) + 1i * randn(size(ftrue)));
	tmp = f.sig * (randn(size(smap)) + 1i * randn(size(smap)));
	zlj = smap .* repmat(ftrue, [1 1 ncoil]) + tmp;
	dB = @(x) 20 * log10(x);
%	f.snr_wrong = dB(norm(ftrue(:)) / norm(yj-ftrue));
%	pr f.snr_wrong
	f.snr = dB(norm(ftrue(:)) / norm(col(yj-ftrue))); % corrected 2012-02-21
	pr f.snr

	slim = [0 1.4];
	im(2, abs(smap), '|smap| true', slim), cbar
	im(3, abs(zlj), '$|$surface coil images $z_{lj}|$'), cbar
	drawnow
end


%% conventional "ratio" estimate of sensitivity map
if ~isvar('shat1'), printm 'ratio'
	for imap = 1:ncoil
		zj = zlj(:,:,imap);
		tmp = zj ./ yj; % usual ratio
		if 1 % set all uncertain map values to median of good ones
			good = abs(zj) > 0.05 * max(abs(zj(:)));
			tmp(~good) = median(abs(tmp(good)));
		end
		shat1(:,:,imap) = tmp;
	end
	clear good zj tmp

	elim = [0 0.25];
	im(4, abs(shat1), '|ratio|', slim), cbar
	im(7, abs(shat1-smap), '|err| ratio', elim), cbar
	alim = [];
	alim = [-pi pi];
	atick = [-3.14 0 3.14];
	im(5, 'hsv', angle(shat1), '$\angle$ ratio', alim), cbar(atick)
	im(8, 'hsv', angle(shat1-smap), '$\angle$ err ratio', alim), cbar(atick)
	drawnow
prompt
end


%% QPLS (quadratically penalized LS) approaches hereafter
if ~isvar('args')
%	f.l2b = -3;
	f.l2b = 2; % log_2(beta) regularization parameter (FWHM = 3.6 pixels)
	args = {zlj, 'order', 2, 'l2b', f.l2b, 'bodycoil', yj};
%	args = {args{:}, 'bodycoil', []}; % test SSoS
%	args = {args{:}, 'bodycoil_default', 'ssos-angle-pca1'};
end


%% cholesky approach that works fine for small 2d images here
if ~isvar('shat2'), printm 'qpls via cholesky'
	cpu etic
	shat2 = mri_sensemap_denoise(args{:}, 'chol', 1, 'niter', 1);
	cpu etoc 'cholesky time'

	im(5, shat2, '|qpls:chol|', slim), cbar
	im(8, abs(shat2-smap), '|err| chol', elim), cbar
	im(6, 'hsv', angle(shat2), '$\angle$ qpls:chol', alim), cbar(atick)
	im(9, 'hsv', angle(shat2-smap), '$\angle$ err chol', alim), cbar(atick)
	drawnow
prompt
end


%% iterative methods hereafter needed for larger problems only, e.g., 3D
%% (but tested here for small 2D problem too)
%% these require remarkably many iterations to approach minimizer...


%% incomplete cholesky method
if ~isvar('shat4'), printm 'qpls via ichol'
	f.niter_ichol = 400;
	cpu etic
	shat4 = mri_sensemap_denoise(args{:}, ...
		'precon', 'ichol', 'niter', f.niter_ichol);
	cpu etoc 'pcg ichol time'
	im(6, abs(shat4), '|qpls:ichol|', slim), cbar
	im(9, abs(shat4-smap), '|err| ichol', elim), cbar
	pr nrms(shat4(:), shat2(:))
prompt
end


%% CG with no preconditioner
if ~isvar('shat3'), printm 'qpls via CG'
	f.niter = 400;
	cpu etic
	shat3 = mri_sensemap_denoise(args{:}, 'niter', f.niter);
	cpu etoc 'pcg time'
	im(6, abs(shat3), '|qpls:iter|', slim), cbar
	im(9, abs(shat3-smap), '|err|', elim), cbar
	pr nrms(shat3(:), shat2(:))
prompt
end


%% examine phase estimates
if 1
	mask = conv2(single(abs(ftrue) > 0), ones(5), 'same') > 0; % dilate
	afun = @(est, tru) angle(est .* conj(tru)) .* mask;
	tfun = @(est, tru) titlef('$\angle$ err, RMS=%.1f$^{\circ}$', ...
		rad2deg(rms(col(masker(afun(est, tru), mask)))));
	im plc 3 4
	im(1, 'hsv', angle(ftrue), '$\angle$ ftrue', alim), cbar(atick)
	im(2, 'hsv', angle(smap), '$\angle$ smap', alim), cbar(atick)
	im(3, 'hsv', angle(zlj), '$\angle$ zlj', alim), cbar(atick)
	im(5, 'hsv', angle(shat1), '$\angle$ ratio', alim), cbar(atick)
	im(9, 'hsv', afun(shat1, smap), '$\angle$ err', alim), cbar(atick)
	tfun(shat1, smap)
	im(6, 'hsv', angle(shat2), '$\angle$ qpls:chol', alim), cbar(atick)
	im(10, 'hsv', afun(shat2, smap), '$\angle$ err', alim), cbar(atick)
	tfun(shat2, smap)
	im(7, 'hsv', angle(shat3), '$\angle$ qpls:iter', alim), cbar(atick)
	im(11, 'hsv', afun(shat3, smap), '$\angle$ err', alim), cbar(atick)
	tfun(shat3, smap)
	im(8, 'hsv', angle(shat4), '$\angle$ qpls:ichol', alim), cbar(atick)
	im(12, 'hsv', afun(shat4, smap), '$\angle$ err', alim), cbar(atick)
	tfun(shat4, smap)
end


return
%% obsolete stuff below here


if 0 && ~isvar('shat2b'), printm 'try precon'
	[shat2a dummy cost2a] = mri_sensemap_denoise(args{:}, ...
		'precon', 1, 'isave', 'all');
	[shat2b dummy cost2b] = mri_sensemap_denoise(args{:}, ...
		'precon', 'diag', 'isave', 'all');
%	plot cost function to see if smaller faster - very small difference?

	clf, plot(0:f.niter, cost2a(:,1), '-o', ...
		0:f.niter, cost2b(:,1), '-x')
	axisx(0,20)
	legend('normal', 'diag')
return
end


%
% old approach (obsolete due to mri_sensemap_denoise)
%
if 0 || ~isvar('shat2'), printm 'regularized method'
	mask1 = true(nx,ny);
	A = diag_sp(yj(mask1));
%	A = diag_sp(ftrue(mask1)); 'testing with ftrue'

%	R = Robject(mask1, 'beta', 2^4);
	R = Robject(mask1, 'beta', 2^f.l2b, 'order', 2, 'distance_power', 2);

	wj = abs(yj);
	wj = median(wj(wj > 0.05 * max(wj(:))));
	W = 1 ./ wj^2; % trick: these weights make the beta "universal"
	if 1 % examine psf
		qpwls_psf(1., R, 1., mask1, 1., 'offset', [8 0]);
		psf = qpwls_psf(A, R, 1., mask1, W, 'offset', [8 0]);
		im(224, psf, 'psf')
	prompt
	end

	for imap = 1:ncoil
		init = shat1(:,:,imap);

		ztmp = zlj(:,:,imap);
		if 0 % test precon (doesn't help!?)
			f.niter = 30;
			M = embed(R.handle_diag(R), mask1);
			M = 1e0 + abs(ftrue).^2 + M;
			% im(M), return
			M = diag_sp(1 ./ M(mask1));
			tmp2 = qpwls_pcg1(init(mask1), A, W, ztmp(:), R.C, ...
				'isave', 'all', 'niter', f.niter, 'precon', M);
			tmp2 = embed(tmp2, mask1);
			tmp1 = qpwls_pcg1(init(mask1), A, W, ztmp(:), R.C, ...
				'isave', 'all', 'niter', f.niter);
			tmp1 = embed(tmp1, mask1);
			cost1 = pwls_cost(tmp1, A, W, ztmp(:), R, mask1);
			cost2 = pwls_cost(tmp2, A, W, ztmp(:), R, mask1);
			plot(0:f.niter, cost1, '-o', 0:f.niter, cost2, '-+')
		return
		end

		tmp = qpwls_pcg1(init(mask1), A, W, ...
			ztmp(:), R.C, 'niter', f.niter);
		shat2(:,:,imap) = embed(tmp, mask1);

		if 0 % cf unif init
			init = ones(nx,ny);
%			ej = zeros(nx,ny); ej(nx/2+1,ny/2+1) = 1;
%			ztmp = ztmp + reshape(A * ej(mask1), size(ztmp));
			tmp = qpwls_pcg1(init(mask1), A, W, ...
				ztmp(:), R.C, 'niter', f.niter);
			shat2u(:,:,imap) = embed(tmp, mask1);
		end
	end
	if 0
		mask0 = conv2(double(ftrue > 0), ones(5), 'same') > 0;
		max_percent_diff(shat2, shat2u)
		mask3 = repmat(mask0, [1 1 ncoil]);
		max_percent_diff(shat2 .* mask3, shat2u .* mask3)
		im(221, shat2u, 'shat2u', slim), cbar
	end
	im(224, shat2, 'Regularized sense map estimates', slim), cbar
prompt
end
%psf = reale(shat2u-shat2, 'warn');
%psf = psf(:,:,1);
%im(223, psf, 'psf'), cbar
%fwhm2(psf)

if 1 % figures for viewing
	im plc 4 1
	im row 1
	im(1, abs(zlj), 'Array coil images', [0 120]), cbar
	im(2, smap, 'True sensivity maps', slim), cbar
	im(3, abs(shat1), 'Ratio sensitivity maps', slim), cbar
	im(4, abs(shat2), 'Regularized sensitivity maps', slim), cbar
return
end

if 1 % figures for publication.  todo: look at real/imag parts
	im clf
	im(smap, 'True sensivity maps', slim), cbar
%	ir_savefig 'fig_mr_sensemap1_map'
	im(abs(zlj), 'Array coil images')
%	ir_savefig 'fig_mr_sensemap1_zl'
	im(abs(shat1), 'Ratio sensitivity maps', slim), cbar
%	ir_savefig 'fig_mr_sensemap1_ratio'
	im(abs(shat2), 'Regularized sensitivity maps', slim), cbar
%	ir_savefig 'fig_mr_sensemap1_reg'
return
end

if 1
	mask3 = 1;
	serr1 = (shat1 - smap) .* mask3;
	serr2 = (shat2 - smap) .* mask3;

	im(121, abs(serr1), 'Error: ratio estimates', elim)%, cbar h
	im(122, abs(serr2), 'Error: regularized estimates', elim)%, cbar h
%	ir_savefig 'fig_mr_sensemap1_err'
return
end
