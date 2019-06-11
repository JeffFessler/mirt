%% mri_example_b0.m
% Example of field-corrected MR reconstruction for a spiral k-space trajectory.
% Copyright 2005-7-29, Jeff Fessler, University of Michigan

%% true zmap
if ~isvar('zmap'), printm 'read fieldmap'
	[fmap, mask_orig] = mri_phase_denoise('test');
	fmap = fmap / (2*pi*2e-3); % 2ms echo time diff.
	[nx, ny] = size(fmap);
	ig = image_geom('nx', nx, 'ny', nx, 'fov', 22); % 22 cm fov
	ig.mask = ellipse_im(ig, [0 0 27*ig.dx 35*ig.dx 0 1], 'oversample', 3) > 0;

	if 1 && im % picture of field map
		t = (conv2(single(ig.mask), ones(3), 'same') >= 1) - ig.mask;
		t = t .* (1 + (-1).^[1:nx]' * (-1).^[1:ny])/2;
		t = t * max(fmap(:));
		if im
			im clf, im(fmap + t)
			title 'Fieldmap and mask outline'
			ir_colorbar_text('Hz')
		end
	end

	zmap = 0 + (2i*pi) * fmap;
prompt
end


%% true image
if ~isvar('xtrue'), printm 'xtrue'
	f.eparam = [
		[0 0 23 29 0 1];
		[0 16 10 10 0 1];
		[0 -16 8 8 0 1];
		[-8 -1 4 8 0 -1];
		[+8 -1 4 8 0 -1];
	];
	f.eparam(:,1:4) = f.eparam(:,1:4) * ig.dx;
	xtrue = ellipse_im(ig, f.eparam, 'oversample', 3);
	xtrue(48:49, 32+[-10:10]) = 2;
	xtrue(16, 32+[-2:2]) = 2;
	xtrue(end/2+1, end/2+1) = 2;

	f.clim = [0 2.5];
	if 1 && im
		t = (conv2(single(ig.mask), ones(3), 'same') >= 1) - ig.mask;
		im clf, im(xtrue + 2.0*t, f.clim), cbar
		title('True image and support outline')
	prompt, clear t
	end
end


%% kspace trajectory
if ~isvar('kspace'), printm 'kspace (slow due to voronoi)'
	N = [nx ny];
	f.traj = 'spiral3';
%	f.traj = 'cartesian';
	f.dens = {'voronoi'}; % todo: use something better!
	[kspace, omega, wi_traj] = mri_trajectory(f.traj, {}, N, ig.fov, f.dens);
	if im
		plot(omega(:,1), omega(:,2), '.')
		axis_pipi, axis square
	prompt
	end
end


%% readout sample times, starting at 0 (single-shot)
if ~isvar('ti')
	f.dt = 5e-6;
	ti = ([1:length(kspace)]-1) * f.dt;
%	ti = ti + .03; warn 'ti shift!' % shift echo for sangwoo test
	f.daq = max(ti) + f.dt;
	printm('readout time: %g ms', 1000 * f.daq)
end



%% "exact" system for generating the data, for less inverse crime
if ~isvar('Ge_zmap'), printm 'Ge_zmap'
	f.basis = 'rect';
%	f.basis = 'dirac*dx';
	Ge_ft = Gmri(kspace, ig.mask, 'exact', 1, 'n_shift', N/2, ...
		'fov', ig.fov, 'basis', {f.basis});
	% trick! adjust DCF wi's to undo the basis effect for CP
	wi_basis = wi_traj ./ Ge_ft.arg.basis.transform;
	if streq(f.basis, 'dirac')
		wi_basis = wi_basis * abs(ig.dx * ig.dy); % trick
	end

	Ge_zmap = feval(Ge_ft.arg.new_zmap, Ge_ft, ti, zmap, {});
	clear yi
end


%% base Gmri object based on NUFFT
if ~isvar('Gn'), printm 'Gn'
	f.nufft = {N, [6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'};
	Gn = Gmri(kspace, ig.mask, 'fov', ig.fov, 'basis', {f.basis}, ...
		'nufft', f.nufft);
end


%% data
if ~isvar('yi'), printm 'yi'
	yt = Ge_zmap * xtrue;
	% add noise
	rng(0)
	yn = randn(size(yt)) + 1i * randn(size(yt));
%	f.snr_db = 50;
	f.snr_db = inf;
	f.scale_noise = norm(yt) / norm(yn) / 10.^(f.snr_db / 20);
	yi = yt + f.scale_noise * yn;
	printm('data rmse = %g, snr = %g', rms(yi-yt), ...
		20*log10(norm(yt(:)) / norm(yi(:)-yt(:))))
	clear xup0 xcp0 yn %xe
end


%% IFFT recon (for cartesian only)
if ~isvar('xifft') && streq(f.traj, 'cartesian'), printm 'xifft'
	xifft = fftshift(ifft2(ifftshift(reshape(yi, ig.dim))));
	if streq(f.basis, 'rect') || streq(f.basis, 'dirac*dx')
		xifft = xifft / abs(ig.dx * ig.dy);
	end
	im plc 1 3
	im(xtrue), cbar
	im(xifft), cbar
	im(xifft-xtrue, 'diff'), cbar
prompt
end


%% "uncorrected" conj. phase recon
if ~isvar('xup0'), printm 'xup0'
	xup0 = Ge_ft' * (wi_basis .* yi);
	xup0 = embed(xup0, ig.mask);
	im clf
	im([xtrue; xup0], 'Uncorrected Conjugate Phase Reconstruction'), cbar
	im plc 1 3
	im(xtrue), cbar
	im(xup0), cbar
	titlef('Uncorrected Conjugate Phase Reconstruction')
	im(xup0-xtrue, 'diff'), cbar
prompt
end


%% slow "exact" conj. phase recon for comparison
if ~isvar('xcp0'), printm 'xcp0'
	xcp0 = Ge_zmap' * (wi_basis .* yi);
	xcp0 = embed(xcp0, ig.mask);
	im plc 1 3
	im(xtrue), cbar
	im(xcp0), cbar
	titlef('Conjugate Phase Reconstruction')
	im(xcp0-xtrue, 'diff'), cbar
prompt
end


%% regularizer
if ~isvar('R'), printm 'R'
	% scale beta by fov^4 since A'A and 2D.
	f.beta = 2^-28 * size(omega,1) * ig.fov^4;
	if streq(f.basis, 'dirac') % needs different beta
		f.beta = f.beta / abs(ig.dx * ig.dy)^2;
	end
	R = Reg1(ig.mask, 'beta', f.beta, 'type_penal', 'mat');

	if 1 % explore resolution
		[psf, var] = qpwls_psf(Ge_ft, R.C, 1, ig.mask);
		im clf, im(psf), cbar
		printm('stddev = %g', sqrt(var * prod(N)))
	prompt
	end
end


%% system model for b0-corrected reconstruction
if ~isvar('Gm'), printm 'Gmri'
	f.L = 6;
	Gm = feval(Gn.arg.new_zmap, Gn, ti, zmap, f.L);
end

if 0 % check approximation accuracy
	max_percent_diff 'Ge_zmap * xtrue(ig.mask)' 'Gm * xtrue(ig.mask)'
return
end


%% CG recon with B0 correction
if ~isvar('xcg1'), printm 'xcg1 iterative'
	f.niter = 15;

	xinit = xcp0;
	xcg1 = qpwls_pcg1(xinit(ig.mask), Gm, 1, yi(:), R.C, 'niter', f.niter);
	xcg1 = embed(xcg1(:,end), ig.mask);
%	im clf, im([xtrue; xcg1], 'xcg1'), cbar
	im plc 1 3
	im(xtrue), cbar
	im(xcg1), cbar
	titlef('CG Reconstruction')
	im(xcg1-xtrue, 'diff'), cbar
prompt
end


%% Toeplitz-based CG
if 0
	if ~isvar('bb'), printm 'bb'
		bb = Gm' * yi(:);
		im clf, im(ig.embed(bb), 'bb'), cbar
		prompt
	end

	if ~isvar('Tm'), printm 'Tm: Gmri gram'
		Tm = build_gram(Gm, 1);
	end

	if ~isvar('xcg2'), printm 'xcg2'
		xcg2 = qpwls_pcg2(xinit(ig.mask), Tm, bb, R.C, 'niter', f.niter);
		xcg2 = ig.embed(xcg2(:,end));
		im clf, im(xcg2, 'xcg2'), cbar
		prompt
	end
end


%% images
if im
	clf
	tmp = [ [xtrue, xcp0]; [xup0, xcg1] ];
	im('notick', abs(tmp), f.clim)
	axis off, title ''
	cbar
	tt = @(x,y,s) ...
		text(x, y, s, 'horizontalalignment', 'center', 'fontsize', 18);
	tt(nx/2, -0.1*ny, 'True');
	tt(3*nx/2, -0.1*ny, sprintf('Uncorrected'));
	tt(nx/2, 2.1*ny, sprintf('Conj. Phase'));
%	tt(nx/2, 2.23*ny, sprintf('L=%d', f.L));
	tt(3*nx/2, 2.1*ny, sprintf('CG-NUFFT'));
	tt(3*nx/2, 2.23*ny, sprintf('L=%d', f.L));
end
