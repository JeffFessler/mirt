%% mri_stack3.m
%|
%| Example illustrating 3D regularized iterative reconstruction for MRI
%| from stack-of-spirals k-space samples.
%|
%| This example does not include field inhomogeneity or relaxation
%| or sensitivity maps.  All of those are possible extensions using
%| existing capabilities of MIRT.
%|
%| Key components below are "kronI" for handling a stack efficiently,
%| and precomputing the iFFT along the kz (stack) direction.
%|
%| For Luis, the number of slices is odd. :)
%|
%| Copyright 2022-08-31, Jeff Fessler, University of Michigan


%% create object
% This is an "honest" simulation where the Fourier data is calculated
% analytically from a continuous-space object model, even though the
% reconstructions are all discrete-space.  So no "inverse crime" here.
if ~isvar('xtrue'), printm 'setup object'

	ig = image_geom('nx', 64, 'nz', 17, ...
		'offsets', 'dsp', ... % (-n/2:n/2-1) for mri
		'fov', [20 20 5.1]); % 20 cm transaxial FOV, 3mm slices
	fov2 = ig.fovs(1:2);
	iy = ceil((ig.ny+1)/2);
	iz = ceil((ig.nz+1)/2);

	xs = mri_objects('fov', ig.fovs, 'test4'); % object model
	% samples of continuous-space, with partial volume in z:
	xtrue = xs.image(ig.xg, ig.yg, ig.zg, 'dx', ig.dx, 'dy', ig.dy, 'dz', ig.dz);

	im plc 3 4
	clim = [0 1] * max(xtrue(:));
	im(xtrue, 'x true', clim), cbar
prompt
end


%% stack-of-spirals k-space trajectory
if ~isvar('kspace'), printm 'trajectory'
%	ktype = 'cartesian';
	ktype = 'spiral3';
	[kspace2, ~, ~] = mri_trajectory(ktype, {}, ig.dim(1:2), fov2); % 2d spiral
	nsample2 = size(kspace2,1); % # of samples in one 2d spiral

	if im % plot trajectory
		omega2 = kspace2 .* (2*pi * [ig.dx ig.dy]);
		im subplot
		plot(omega2(:,1), omega2(:,2), '.')
		titlef('One spiral with %d samples', size(omega2,1))
		axis_pipi, axis square
	end

	kz = ig.kz; % DFT sampling in kz
	kspace      = col(repmat(kspace2(:,1), ig.nz, 1));
	kspace(:,2) = col(repmat(kspace2(:,2), ig.nz, 1));
	kspace(:,3) = col(repmat(kz, [nsample2 1]));

	if 0 % check
		% plot3(kspace(:,1), kspace(:,2), kspace(:,3), '-.')
		[kspace3, omega3, wi_traj3] = mri_trajectory('spiral3', {}, ig.dim, ig.fovs, f.dens);
		clf, plot([kspace(:,3) kspace3(:,3)])
		max_percent_diff(kspace3, kspace), return
	end
prompt
end


%% simulate k-space data
if ~isvar('yi_stack'), printm 'k-space data'
	printm 'setup data, based on continuous-space object model'
	ytrue = xs.kspace(kspace(:,1), kspace(:,2), kspace(:,3));
	ytrue = reshape(ytrue, nsample2, ig.nz); % remind us it is a stack!

	% add noise, if desired
	rng(0)
	yi_stack = ytrue + 0 * (randn(size(ytrue)) + 1i * randn(size(ytrue)));
	pr snr(ytrue, yi_stack - ytrue)
end


%% Take IFFT in kz once and for all
if ~isvar('yi_fftz')
	% note: "2" is kz! and we need 1/dz for proper scaling
	fftz = @(y) fftshift(ifft(ifftshift(y,2), [], 2), 2) / ig.dz;
	yi_fftz = fftz(yi_stack);
end


%% create Gmri object
if ~isvar('A3'), printm 'system model'
	N = ig.dim(1:2);
	nufft_args = {N, 6*ones(size(N)), 2*N, N/2, 'table', 2^10, 'minmax:kb'};
	mask2 = true(N);
	clear N
	f.basis = {'rect'};
%	f.basis = {'dirac*dx'};
	A2 = Gmri(kspace2, mask2, ...
		'fov', fov2, 'basis', f.basis, 'nufft', nufft_args);
	A3 = kronI(ig.nz, A2); % system model for stack!
end


%% Density compensation factors (DCF)
if ~isvar('wi_traj2'), printm 'DCF'

	% voronoi
	wi_voro2 = ir_mri_density_comp(kspace2, 'voronoi', 'fix_edge', 0);
	wi_max = 1.05 / prod(fov2);
	wi_voro2 = min(wi_voro2, wi_max);

	% pipe & menon
	wi_pipe2 = ir_mri_density_comp(kspace2, 'pipe', 'G', A2.arg.Gnufft, ...
		'arg_pipe', {'fov', fov2, 'niter', 60});

	if im % plot DCF
		im subplot
		plot([1 nsample2], wi_max * [1 1], 'm-', ...
			1:nsample2, wi_voro2, '.', ...
			1:nsample2, wi_pipe2, '.')
		titlef('DCF for one 2D spiral')
		xlim([1 nsample2]), xtick([1 nsample2])
		legend('wi max', 'DCF voronoi', 'DCF Pipe', ...
			'location', 'southeast'), drawnow
	end

	% trick! undo basis effect and pixel size from FT of sinc3(x/d)
	basis_correction = (ig.dx * ig.dy) ./ A2.arg.basis.transform;

	if im % check 2D PSF
		psf_voro2 = embed(A2' * (basis_correction .* wi_voro2), mask2);
		cri = @(x) cat(3, real(x), imag(x));
		im(cri(psf_voro2))
		xlabelf('Re(sum) is %.2f', real(sum(psf_voro2(:))))

		psf_pipe2 = embed(A2' * (basis_correction .* wi_pipe2), mask2);
		im(cri(psf_pipe2))
		xlabelf('Re(sum) is %.2f', real(sum(psf_pipe2(:))))

		im subplot
		tmp = [psf_voro2(:,iy) psf_pipe2(:,iy)];
		plot(ig.x, real(tmp), 'b.-', ig.x, imag(tmp), 'r.-')
		legend('real voronoi', 'real pipe', 'imag voronoi', 'image pipe')
	end

	% use Pipe&Menon because it has better-scaled sum
	wi_traj2 = wi_pipe2; clear wi_pipe2 wi_voro2
prompt
end


%% Examine PSF of conjugate-phase reconstruction using DCF
if ~isvar('psf3')
	wi_basis2 = wi_traj2 ./ A2.arg.basis.transform; % trick! undo basis effect
	fake_data = wi_basis2 .* ones(size(yi_fftz)); % broadcast!
	fake_data = fake_data * ig.dx * ig.dy * ig.dz; % trick from FT of sinc3(x/d)
	fake_data = fftz(fake_data);
	psf3 = A3' * fake_data; % broadcast!
%	psf3 = embed(psf3, ig.mask);
	im('mid3', abs(psf3), '3D PSF slices'), cbar

	if im
		im subplot
		plot(ig.x, real(psf3(:,iy,iz)), 'b.-', ig.x, imag(psf3(:,iy,iz)), 'r.-')
		legend('real', 'imag')
	end
prompt
end


%% Conjugate phase reconstruction
if ~isvar('xcp'), printm 'conj. phase reconstruction'
	minmax(A2.arg.basis.transform)
	wi_basis2 = wi_traj2 ./ A2.arg.basis.transform; % trick! undo basis effect
	minmax(wi_basis2)

	xcp = A3' * (wi_basis2 .* yi_fftz); % slice-wise broadcast
%	xcp = embed(xcp, mask);
	im(abs(xcp), '|Conj. Phase| Recon'), cbar
	xlabelf('NRMSE = %.1f%%', 100*nrms(xcp, xtrue))

	if im
		im subplot
		plot(	ig.x, squeeze(xtrue(:,iy,iz)), '.-', ...
			ig.x, squeeze(real(xcp(:,iy,iz))), '-o')
		ir_legend({'\x true', 'Re(CP)'}), titlef('Profiles')
		xlabelf('x [cm]')
	end
prompt
end


%% regularizer
if ~isvar('R'), printm 'regularizer'
%	beta = 2^-7 * size(kspace,1); % good for quadratic 'dirac'
	beta = 2^-21 * size(kspace,1); % good for quadratic 'rect'
	R = Reg1(ig.mask, 'beta', beta * [1 1 0]); % trick: none in z!
	if 1 % check resolution: [1.14 1.14 1] for 'rect'
		qpwls_psf(A3, R, 1, ig.mask, 1, 'fwhmtype', 'profile');
%	return
	end
end


%% PCG
if ~isvar('xpcg'), printm 'PCG with quadratic regularizer'
	niter = 10;
	xpcg = qpwls_pcg1(xcp(ig.mask), A3, 1, yi_fftz(:), R.C, 'niter', niter);
	xpcg = ig.embed(xpcg);
	im(abs(xpcg), '|PCG quad|', clim), cbar
	xlabelf('NRMSE = %.1f%%', 100*nrms(xpcg, xtrue))

	if im
		im subplot
		plot(	ig.x, squeeze(xtrue(:,iy,iz)), '.-', ...
			ig.x, squeeze(real(xpcg(:,iy,iz))), '-+')
		ir_legend({'\x true', 'Re(PCG)'}), titlef('Profiles')
		xlabelf('x [cm]')
	end
end
