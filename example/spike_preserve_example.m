%| spike_preserve_example.m
%| Examine how much edge-preserving regularization preserves spikes,
%| in denoising application:
%| min_x 1/2 |y - x|^2 + \beta R(x)
%|
%| Copyright 2010-08-05, Jeff Fessler and Won Huh, University of Michigan

if ~isvar('yi')
	nx = 2^7;
	ny = nx;
	xtrue = zeros(nx, ny);
	iy = 1:ny;
	ramp = 1 + 150*iy/ny;
	xtrue(nx/4,iy) = 1 + 100*iy/ny;
	xtrue(nx/2:3*nx/4,iy) = 20;
	im plc 2 2
	clim = [-5 105];
	im(1, xtrue, clim)

	% noisy data
	rng(30)
	sig = 8;
	yi = xtrue + sig * randn(nx,ny);
	im(2, yi, clim)

	mask = true(nx,ny);
end

	f.delta = 2^-1;
	R = Reg1(mask, 'type_denom', 'matlab', ...
		'offsets', [1], ... % trick for 1D only
		'beta', 2^7, 'pot_arg', {'hyper3', f.delta});

	A = 1; % identity (denoising)
	f.niter = 200;
	xh = pwls_pcg1(yi(:), A, 1, yi(:), R, 'niter', f.niter);
	xh = embed(xh, mask);
	im(3, xh, clim)

	im subplot 4
%	plot([xh(:,[nx/2 end]) xtrue(:,[nx/2 end])])
%	plot(iy, ramp, '-', iy, xh(nx/4,:), '.')
	plot(ramp, xh(nx/4,:), '.', ramp, ramp, '--')
	xlabel 'true spike intensity'
	ylabel 'estimated spike intensity'
	axis tight, axis square
	title 'weak spikes are suppressed more, cf soft thresholding'

printm 'todo: open challenge is to find where knee in curve is'
printm 'as a function of sigma, delta, beta'
