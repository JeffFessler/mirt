%% mri_example3.m
%|
%| Example illustrating edge-preserving iterative reconstruction for MRI.
%| Specifically, show how it can "recover" some high spatial frequencies.
%|
%| Copyright 2010-07-22, Jeff Fessler, University of Michigan

%% functions for true object and its (analytical) spectrum - no inverse crime!
if ~isvar('xtrue'), printm 'setup object'
	f.fov = 256; % FOV in mm

	f.nxc = 2^6; % coarse sampling
	f.nxf = 2^8; % fine sampling
	igc = image_geom('nx', f.nxc, 'fov', f.fov);
	igf = image_geom('nx', f.nxf, 'fov', f.fov);

	obj = mri_objects('circ2', [0 0 f.fov/3 1]);
	xtrue = obj.image(igf.xg, igf.yg);

	im plc 4 2
	clim = [-0.1 1.2];
	im(1, igf.x, igf.y, xtrue, 'x true (fine)', clim), cbar
	xlabel 'x [mm]', ylabel 'y [mm]'

	Xf = obj.kspace(igf.ug, igf.vg);

	iv = find(igf.v == 0);
	ff = @(f) log(abs(f));

	f.af = [-0.5 0.5 -2 11];
	im subplot 2
	plot(igf.u, ff(Xf(:,iv)))
	xlabelf('\kx [1/mm]')
	ylabelf('$\log(|X(\kx,0)|$')
	titlef 'True spectrum'
	axis(f.af)
prompt
end


%% standard iFFT recon
if ~isvar('xfft'), printm 'xfft'
	Xc = obj.kspace(igc.ug, igc.vg);
	ic = find(igc.v == 0);
	im(3, Xc)

	im subplot 4
	plot(igc.u, ff(Xc(:,ic)))
	xlabel 'u [1/mm]', ylabel 'log(|X(u,0)|)'
	title 'x spectrum - kspace samples'
	axis(f.af)

	tmp = igf.zeros;
	tmp(1:igc.nx,1:igc.ny) = Xc;
	tmp = circshift(tmp, -[igc.nx igc.ny]/2);
	im(3, tmp)

	xfft = fftshift(ifft2(tmp)) / abs(igf.dx * igf.dy);
	xfft = reale(xfft, 0.02);
	im(3, igf.x, igf.y, xfft, clim), cbar
	xlabel 'x [mm]', ylabel 'y [mm]'
	title 'iFFT recon'
prompt
end


%% system matrix
if ~isvar('A'), printm 'A'
%	igf.mask = igf.circ > 0;
	im(8, igf.mask)

	samp = igf.zeros;
	samp(1:igc.nx,1:igc.ny) = 1;
	samp = circshift(samp, -[igc.nx igc.ny]/2);
	samp = logical(samp);
	im(8, samp)

	A = Gdft('mask', igf.mask, 'samp', samp, 'ifftshift', true);

	yb = fftshift(Xc);
	yb = yb(:);

	if 1 % check it: todo 1/2 pixel shift
		duv = 1 / abs(igf.nx * igf.dx * igf.ny * igf.dy);
		tmp = A' * yb * duv; % conj. phase reconstruction
		xcp = igf.embed(tmp);
		xcp = reale(xcp, 0.02);
		im(8, igf.x, igf.y, xcp, clim), cbar
	end
prompt
end


%% CG quadratic regularizer
if ~isvar('xq'), printm 'CG with quadratic regularizer'
	beta = 2^-2;
	Rq = Reg1(igf.mask, 'beta', beta);

	if 0 % example PSF
		qpwls_psf(A, R, 1, igf.mask);
	return
	end

	f.niter = 10;
	xq = qpwls_pcg1(xcp(igf.mask), A, 1, yb(:), Rq.C, ... 
		'niter', f.niter);
	xq = igf.embed(xq(:,end));
	xq = reale(xq, 'warn');

	im(5, igf.x, igf.y, xq, 'CG quad', clim), cbar
	xlabel 'x [mm]', ylabel 'y [mm]'

	Xq = fftshift(fft2(xq));
	im subplot 6
	plot(igf.u, ff(Xq(:,iv)))
	axis(f.af)
	xlabelf('\kx [1/mm]'), ylabelf('$\log(|X(\kx,0)|$')
	title 'CG quad spectrum'
prompt
end


%% CG edge preserving
if ~isvar('xh'), printm 'CG with edge-preserving regularizer'
	R = Reg1(igf.mask, 'beta', 2^3*beta, 'pot_arg', {'hyper3', 0.1}, ...
			'type_penal', 'mat');
	xh = pwls_pcg1(xcp(igf.mask), A, 1, yb(:), R, ...
		'niter', 1*f.niter);
	xh = igf.embed(xh);
	xh = reale(xh, 'warn');

	im(7, igf.x, igf.y, xh, 'CG edge', clim), cbar
	xlabel 'x [mm]', ylabel 'y [mm]'

	Xe = fftshift(fft2(xh));
	im subplot 8
	plot(igf.u, ff(Xe(:,iv)))
	axis(f.af)
	xlabelf('\kx [1/mm]'), ylabelf('$\log(|X(\kx,0)|$')
	title 'CG edge spectrum'
end


if 0
	clf
	%im subplot 4
	plot([xfft(:,end/2) xq(:,end/2), xh(:,end/2)])
	legend('fft', 'quad', 'edge')
end
