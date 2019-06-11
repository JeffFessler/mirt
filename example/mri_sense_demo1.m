%% mri_sense_demo1.m
%|
%| Example illustrating regularized iterative reconstruction for parallel MRI
%| (sensitivity encoding imaging or SENSE), from nonuniform k-space samples.
%| (This example does not include field inhomogeneity or relaxation.)
%|
%| Copyright 2006-4-18, Jeff Fessler, University of Michigan

%% make coil sensitivity maps
if ~isvar('smap'), printm 'sense maps'
	ig = image_geom('nx', 64, 'fov', 250); % 250 mm FOV
	ig.mask = ig.circ(ig.dx * (ig.nx/2-2), ig.dy * (ig.ny/2-1)) > 0;
	f.ncoil = 4;
	smap = mri_sensemap_sim('nx', ig.nx, 'ny', ig.ny, 'dx', ig.dx, ...
		'rcoil', 120, 'ncoil', f.ncoil, 'chat', 1);
prompt
end


%% true object (discrete because of discrete sense map)
if ~isvar('xtrue'), printm 'true object'
	xtrue = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
	clim = [0 8];
	im plc 2 2
	im(1, ig.x, ig.y, xtrue, '\x true', clim), cbar
	im(4, ig.mask, 'mask')
prompt
end


%% trajectory and system matrix
if ~isvar('Am'), printm 'system matrix A'
	f.traj = 'spiral1';
	f.dens = 'voronoi';

	N = [ig.nx ig.ny];
	[kspace, omega, wi_traj] = mri_trajectory(f.traj, {}, ...
		N, ig.fov, {f.dens});

	% create Gnufft class object
	J = [6 6];
	nufft_args = {N, J, 2*N, N/2, 'table', 2^10, 'minmax:kb'};
	Am = Gmri(kspace, ig.mask, ...
		'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);

	if im
		im subplot 2
		plot(omega(1:5:end,1), omega(1:5:end,2), '.')
		axis_pipi, axis square
		titlef('%s: %d', f.traj, size(omega,1))
	end
prompt
end

wi_basis = wi_traj ./ Am.arg.basis.transform;


%% SENSE system object
if ~isvar('Ab'), printm 'Ab object with sense maps within'
	for ic=1:f.ncoil
		tmp = smap(:,:,ic);
		tmp = Gdiag(tmp(ig.mask), 'mask', ig.mask);
		Ac{ic} = Am * tmp; % cascade
	end
	Ab = block_fatrix(Ac, 'type', 'col'); % [A1; A2; ... ]
end


%% noisy data
if ~isvar('yi'), printm 'data yi'
	ytrue = Ab * xtrue(ig.mask);

	% add noise
	rng(0)
	yi = ytrue + 0 * randn(size(ytrue));
	% todo: visualize data...
end


%% CP recon
if ~isvar('xcp'), printm 'conj. phase reconstruction'
	xcp = zeros(ig.nx, ig.ny, f.ncoil, 'single');
	y4 = reshape(yi, [], f.ncoil);
	for ic=1:f.ncoil
		tmp = Am' * (wi_basis .* y4(:,ic));
		xcp(:,:,ic) = ig.embed(tmp);
	end

	im(2, abs(xcp), 'Conj. Phase Recon for each coil'), cbar
prompt
end


%% SSoS
if ~isvar('xssos'), printm 'sqrt sum-of-squares reconstruction'
	xssos = sqrt(sum(abs(xcp).^2, 3));
	xssos = ir_wls_init_scale(1, xtrue, xssos); % cheat trick for scaling
	im(3, abs(xssos), 'sum-of-squares recon'), cbar
	xlabelf('NRMSE %.1f\%%', 100*nrms(xssos(:), xtrue(:)))
prompt
end


%% regularizer
if ~isvar('R'), printm 'regularizer'
	f.beta = 2^11;
	R = Reg1(ig.mask, 'beta', f.beta, 'type_penal', 'mat'); % complex
	if 1
		psf = qpwls_psf(Ab, R, 1, ig.mask, 1, 'offset', [0 0]);
		im(4, psf)
	end
prompt
end


%% PCG
if ~isvar('xpcg'), printm 'PCG with quadratic penalty'
	f.niter = 10;
	xpcg = qpwls_pcg(xssos(ig.mask), Ab, 1, yi(:), 0, R.C, 1, f.niter);
	xpcg = ig.embed(xpcg(:,end)); % convert last vector to image for display

	im(4, abs(xpcg), '|\x| PCG quad', clim), cbar
	xlabelf('NRMSE %.1f\%%', 100*nrms(xpcg(:), xtrue(:)))
end
