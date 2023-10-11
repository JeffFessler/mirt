%% mri_example_3d.m
%|
%| Example illustrating 3D regularized iterative reconstruction for MRI
%| from nonuniform k-space samples.
%|
%| This example does not include field inhomogeneity or relaxation
%| or sensitivity maps.  All of those are possible extensions using
%| existing capabilities of MIRT.
%|
% Copyright 2004-12-16, Jeff Fessler, University of Michigan

%% create object
% This is an "honest" simulation where the Fourier data is calculated
% analytically from a continuous-space object model, even though the
% reconstructions are all discrete-space.  So no "inverse crime" here.
if ~isvar('xtrue'), printm 'setup object'

	nx = 64; nz = 32; % original
	ig = image_geom('nx', nx, 'nz', nz, ...
		'offsets', 'dsp', ... % (-n/2:n/2-1) for mri
		'fov', [20 20 6]); % 20 cm transaxial FOV
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


%% k-space trajectory
if ~isvar('kspace'), printm 'trajectory'
% todo: 3d radial (kushball)
%	f.traj = 'radial'; % stack-of-stars
%	f.traj = 'cartesian';
	f.traj = 'spiral3'; % stack-of-spirals
	[kspace, ~, wi_simple] = mri_trajectory(f.traj, {}, ig.dim, ig.fovs, {});

	if im % plot trajectory
		omega = kspace .* (2*pi * ig.deltas);
		im subplot
		plot3(omega(:,1), omega(:,2), omega(:,3), '.')
		titlef('%s with %d samples', f.traj, size(omega,1))
		axis_pipi, axis square
	end
prompt
end


%% simulate k-space data
if ~isvar('yi'), printm 'k-space data'
	printm 'setup data, based on continuous-space object model'
	ytrue = xs.kspace(kspace(:,1), kspace(:,2), kspace(:,3));

	% add noise, if desired
	rng(0)
	yi = ytrue + 0 * (randn(size(ytrue)) + 1i * randn(size(ytrue)));
	pr snr(ytrue, yi - ytrue)

	dbfun = @(y) max(20*log10(abs(y) / max(abs(y(:)))), -80);
	if streq(f.traj, 'cartesian')
		im(reshape(dbfun(yi), ig.dim)), cbar
		titlef('$\log |y_i|$')
	end
end


%% sanity check in cartesian case
if 0
	ytmp = fftshift(fftn(fftshift(xtrue))); % cartesian
	ytmp = ytmp * abs(ig.dx * ig.dy * ig.dz);
	im(abs(ytmp), 'FFT(xtrue)'), cbar

	y3 = reshape(yi, ig.dim);
	im(abs(y3-ytmp), 'fft err'), cbar
	max_percent_diff(y3, ytmp)

	xtmp = fftshift(ifftn(fftshift(y3)));
	xtmp = xtmp / abs(ig.dx * ig.dy * ig.dz);
	im(xtmp), cbar
	max_percent_diff(xtrue, xtmp)

	im subplot
	plot(	ig.z, squeeze(xtrue(end/2+1,end/2+1,:)), '-o', ...
		ig.z, squeeze(xtmp(end/2+1,end/2+1,:)), '-x')
	ir_legend({'\x true', 'iFFT'})
return
end


%% create Gmri object
if ~isvar('A3'), printm 'system model'
	N = ig.dim;
	nufft_args = {N, 6*ones(size(N)), 2*N, N/2, 'table', 2^10, 'minmax:kb'};
	mask = true(N);
	clear N
	f.basis = {'rect'};
%	f.basis = {'dirac*dx'};
	A3 = Gmri(kspace, mask, ...
		'fov', ig.fovs, 'basis', f.basis, 'nufft', nufft_args);
end


%% Density compensation factors (DCF)
if ~isvar('wi_traj'), printm 'DCF'
	if 1
		f.dens = 'simple';
		wi_traj = wi_simple; % fast
	else
%		f.dens = 'voronoi';
		f.dens = 'pipe'; % pipe & menon
		wi_traj = ir_mri_density_comp(kspace, f.dens, 'G', A3.arg.Gnufft, ...
			'arg_pipe', {'fov', ig.fovs, 'niter', 60});
	end
	wi_max = 1.05 / prod(ig.fovs);
	if 1, minmax(wi_traj), end

	if im % plot DCF
		nsample = size(kspace,1);
		im subplot
		plot([1 nsample], wi_max * [1 1], 'm-', ...
			1:nsample, wi_traj, '.');
		titlef('DCF')
		xlim([1 nsample]), xtick([1 nsample])
		legend('wi max', ['DCF ' f.dens], ...
			'location', 'southeast'), drawnow
	end

	% trick! undo basis effect and pixel size from FT of sinc3(x/d)
	basis_correction = prod(ig.deltas) ./ A3.arg.basis.transform;

	% verify that DCF is identical for all spirals in a stack-of-spirals
	if im && streq(f.traj, 'spiral3')
		im subplot
		tmp = reshape(wi_traj, [], ig.nz);
		plot(tmp, '.'), titlef('DCFs for stack')
	end
prompt
end


%% Examine PSF of conjugate-phase reconstruction using DCF
if ~isvar('psf3')
	wi_basis = wi_traj ./ A3.arg.basis.transform; % trick! undo basis effect
	psf3 = A3' * wi_basis;
	psf3 = embed(psf3, mask);
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
	minmax(A3.arg.basis.transform)
	wi_basis = wi_traj ./ A3.arg.basis.transform; % trick! undo basis effect
	minmax(wi_basis)

	xcp = A3' * (wi_basis .* yi);
	xcp = embed(xcp, mask);
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
	R = Reg1(ig.mask, 'beta', beta);
	if 0 % check resolution: [1.17 1.17 1.01] for 'dirac'
		psfq = qpwls_psf(A3, R, 1, ig.mask, 1, 'fwhmtype', 'profile');
		if 0 % for haowei
			ii = 1 + (-7:7);
			im clf; jim(psfq(end/2+ii, end/2+ii, end/2+ii))
			psf2 = embed(A3' * (A3 * ig.unitv), ig.mask);
			im clf; jim(psf2(end/2+ii, end/2+ii, end/2+ii))
			Psf2 = fftshift(fftn(fftshift(psf2)));
			plot(log10(abs(Psf2(:,end/2+1,end/2+1))))
		end
	return
	end
end


%% PCG
if ~isvar('xpcg'), printm 'PCG with quadratic regularizer'
	niter = 10;
	xpcg = qpwls_pcg1(xcp(ig.mask), A3, 1, yi(:), R.C, 'niter', niter);
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
