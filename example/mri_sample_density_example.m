%% mri_sample_density_example.m
%|
%| (work in progress)
%|
%| Evaluate various sample density compensators
%|
%| Copyright 2003-7-29, Jeff Fessler, University of Michigan


%% create k-space samples, in units of 1/mm
if ~isvar('kspace'), printm 'setup kspace'
	fov = 256; % [mm] typical brain FOV
	N0 = 64; % nominal image size
	kmax = N0/2*(1/fov); % for display axes

%	ktype = 'cartesian';
%	ktype = 'random';
	ktype = 'spiral';

	switch ktype
	case 'random'
		rng(0)
		kspace = (rand(N0*N0, 2)-0.5)*N0/fov;	% random k-space
	case 'spiral'
		t = linspace(0, N0/2*2*pi, N0^2)';	% crude spiral:
		kspace = N0/2*(1/fov)*[cos(t) sin(t)] .* (t(:,[1 1]) / max(t));
	case 'cartesian'
		k1 = [-N0/2:N0/2-1]/fov;		% cartesian
		[kk1, kk2] = ndgrid(k1, k1);
		kspace = [kk1(:), kk2(:)];
		1/diff(kspace(1:2,1))
	end, clear t k1 kk1 kk2

	if im
		im plc 3 3
		im subplot 1
		plot(kspace(:,1), kspace(:,2), '.')
		axis(1.1*[-1 1 -1 1]*kmax), axis square
		xlabelf '\kx [1/mm]', ylabelf '\ky [1/mm]'
		titlef('%d %s k-space samples', size(kspace,1), ktype)
	end
prompt
end


%% true object and analytical k-space data
if ~isvar('xtrue'), printm 'setup object'
	% display images with many pixels...
	x1d = [-N0/2:N0/2-1] / N0 * fov;
	[x1dd, x2dd] = ndgrid(x1d, x1d);

	% parameter units all in [mm]
	obj = mri_objects('case1');
	xtrue = obj.image(x1dd, x2dd);
	ytrue = obj.kspace(kspace(:,1), kspace(:,2));
	clear x1dd x2dd obj

	clim = [0 2];
	im(2, x1d, x1d, xtrue, 'x true', clim), cbar

	% add noise
	rng(0)
	yd = ytrue + 0 * randn(size(ytrue));

	% lazy kspace data gridding for display
	[xg, yd_g, x1g, k1g] = mri_grid_linear(kspace, yd, N0, fov);

	pr 'imax(yd_g,2)'
	im(3, k1g{1}, k1g{2}, abs(yd_g), '$|y_i|$'), cbar
	im(4, x1g{1}, x1g{2}, abs(xg), '|\x| "gridding"', clim), cbar

	mask = true(N0);
prompt
end, clear yd_g k1g


%% create Gnufft class object
if ~isvar('A'), printm 'setup Gnufft object'
	omega = 2*pi*kspace*fov/N0;
%	disp(minmax(omega(:,1)))
	A = Gnufft({omega, [N0 N0], [6 6], 2*[N0 N0], [N0/2 N0/2]});
%	A = Gnufft({omega, [N0 N0], [8 8], 4*[N0 N0], [N0/2 N0/2], 'kaiser'});
prompt
end


%% gridding (without DCF)
if ~isvar('xgr.unif')
	xgr.unif = embed(A' * yd, mask);
	im(5, x1g{1}, x1g{2}, abs(xgr.unif), 'no comp'), cbar
prompt
end


%% Pipe and Menon style DCF
if ~isvar('wt.pipe'), printm 'pipe and menon'
	w = ones(size(yd));
	P = A.arg.st.p;
	for ii=1:20
		tmp = P * (P' * w);
		w = w ./ real(tmp);
	end
	minmax(tmp)
	wt.pipe = w;
	if im
		im subplot 7
		plot(wt.pipe)
		title 'DCF'
	end
end


%% gridding with Pipe&Menon DCF
if ~isvar('xgr.pipe'), printm 'pipe&menon gridding'
	scale = A.arg.st.sn(end/2,end/2)^(-2) / fov^2 / prod(A.arg.st.Kd) * N0^2;
	w = wt.pipe * scale;
	xgr.pipe = embed(A' * (w .* yd), mask);
	if im
		im(8, x1g{1}, x1g{2}, real(xgr.pipe), 'Pipe+Menon', clim), cbar
		xlabelf('NRMSE %.1f\%%', 100*nrms(xgr.pipe(:), xtrue(:)))
	end
	pr 'sum(xgr.pipe) / sum(xtrue)'
end


%% reconstruct by PCG
if ~isvar('xpcg'), printm 'PCG with quadratic penalty'

	niter = 20;
	beta = 2^-10 * size(omega,1);

	R = Reg1(mask, 'beta', beta, 'type_penal', 'mat');
	ytmp = yd(:) * (1/fov)^2 * N0^2; % scaling! todo: why?
	xiter = qpwls_pcg(xg(:), A, 1, ytmp(:), 0, R.C, 1, niter);
	xpcg = embed(xiter(:,end), mask);
	im(3, x1g{1}, x1g{2}, abs(xpcg), '|\x| pcg quad'), cbar
%prompt
end

% ir_savefig fig_mri_sample_density
