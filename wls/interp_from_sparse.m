% interp_from_sparse.m
% Using PWLS algorithms to "interpolate" from sparsely sampled data.
% Work in progress - results are mysteriously bad.  todo: look in 1D?
% Copyright 2010-1-7, Jeff Fessler, University of Michigan

if ~isvar('yi'), printm 'data'
	ig = image_geom('nx', 110, 'ny', 128, 'dx', 1);
	xtrue = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
	im plc 2 3
	im(1, xtrue, 'xtrue (high res)')

	% generate sparse "down sampling" sampling pattern
	samp = false(ig.nx, ig.ny);
	down = 2;
	ix = floor((down+1)/2):down:ig.nx;
	iy = floor((down+1)/2):down:ig.ny;
	samp(ix,iy) = true;

	yi = samp .* xtrue;
	im(2, yi, 'yi (sampled)')
prompt
end


if ~isvar('xc'), printm 'basic interpolation'
	x0 = interpn(ix, iy, yi(ix,iy), [1:ig.nx]', 1:ig.ny, 'nearest');
	x0(isnan(x0)) = 0;
	im(3, x0, 'x0, 0-order hold')

	xc = interpn(ix, iy, yi(ix,iy), [1:ig.nx]', 1:ig.ny, 'cubic');
	xc(isnan(xc)) = 0;
	im(4, xc, 'xc, cubic interpolation')
prompt
end


if ~isvar('Rn'), printm 'R, regularizer'
	Rq = Robject(ig.mask, 'type_denom', 'matlab', ...
		'potential', 'quad', 'beta', 2^(-6));
	f.pot = 'huber';
%	f.pot = 'hyper3';
	Rn = Reg1(ig.mask, 'type_denom', 'matlab', ...
		'pot_arg', {f.pot, 0.3}, 'beta', 2^(-7));
prompt
end


if ~isvar('xq'), printm 'xq: quadratic regularization'
	A = Gdiag(single(samp));
	W = Gdiag(single(samp));
	f.niter = 400;
	xinit = x0;
	xq = qpwls_pcg1(xinit(ig.mask), A, W, yi(:), Rq.C, 'niter', f.niter);
	xq = ig.embed(xq);
	im(5, xq, 'QPWLS')
prompt
end


if 0 % check cost: indeed xq is lowest
	cost = @(x) pwls_cost(x(ig.mask), A, W, yi(:), Rq);
	cost(xtrue);
	cost(x0);
	cost(xq);
return
end


if ~isvar('xs'), printm 'xs: nonquadratic edge-preserving penalty'
	f.niter2 = 150;
	xs = pwls_pcg1(xinit(ig.mask), A, W, yi(:), Rn, 'niter', f.niter2);
	xs = ig.embed(xs);
	im(6, xs, 'PWLS')
prompt
end

if 1
	cost = @(x) pwls_cost(x(ig.mask), A, W, yi(:), Rn);
	cost(xtrue);
	cost(x0);
	cost(xs);
end

if im
	rms = @(x) sqrt(mean(col(abs(x-xtrue)).^2));
	f.r0 = rms(x0);
	f.rc = rms(xc);
	f.rq = rms(xq(:,:,end));
	f.rs = rms(xs(:,:,end));
	pr f

	im pl 2 3
	clim = [0 6];
	im(1, xtrue, clim, 'truth'), cbar
	im(2, yi, clim, 'yi'), cbar
	axis([1 ig.nx 1 ig.ny])
	im(3, x0, clim, 'x0, zero-order'), cbar
	xlabel(num2str(f.r0))
	im(4, xc, clim, 'xc, cubic'), cbar
	xlabel(num2str(f.rc))
	im(5, xq(:,:,end), clim, 'QPWLS'), cbar
	xlabel(num2str(f.rq))
	im(6, xs(:,:,end), clim, 'PWLS edge preserving'), cbar
	xlabel(num2str(f.rs))
end
