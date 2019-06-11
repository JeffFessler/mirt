% ir_example_ct_lir1
% example of local impulse response (LIR) in CT
% Copyright 2012-07-28, Jeff Fessler, University of Michigan

if ~isvar('A'), printm 'setup geometry, image, sinogram'
	f.down = 4;
	ig = image_geom('nx', 512, 'fov', 30, 'down', f.down);
	ig.mask = ig.circ > 0;
	sg = sino_geom('par', 'nb', ig.nx, 'na', ig.nx, ...
			'dr', ig.dx, 'strip_width', 'd');
	sg.plot(ig);

	ell = [0 0 10 10 0 0.2;
                 0  6 3 3 0  0.2;
                 6  0 3 3 0  0.1
                 0 -6 3 3 0 -0.1;
                -6  0 3 3 0 -0.2;
		];
        xtrue = ellipse_im(ig, ell, 'oversample', 2);

%	A = Gtomo2_dscmex(sg, ig);
	A = Gtomo2_wtmex(sg, ig, 'nthread', jf('ncore'));
	ytrue = A * xtrue;
	wi = exp(-ytrue); % ideal Poisson weighting

	im plc 2 3
	clim = [0 0.4];
	im(1, xtrue, 'x', clim), cbar
	im(2, ytrue, 'y: ideal sinogram'), cbar
	im(3, wi, 'w: ideal weighting'), cbar
%	ir_savefig ir_example_ct_lir1_x_y_w
prompt
end


if ~isvar('yp'), printm 'yp'
	ix = round(ell(2:end, 1) / ig.dx + (ig.nx-1)/2);
	iy = round(ell(2:end, 2) / ig.dy + (ig.ny-1)/2);
	xp = ig.zeros;
	for ii=1:numel(ix)
		xp = xp + ig.unitv(ix(ii), iy(ii));
	end
	f.amp = 0.1;
	xp = f.amp * xp ;
	im(2, xtrue + xp)

	yp = A * xp;
	climp = [0 f.amp];
	im(3, xp, climp), cbar
	im(4, yp)
end

if ~isvar('fbp'), printm 'fbp'
	tmp = fbp2(sg, ig);
	fbp = fbp2(yp, tmp, 'window', 'boxcar,0.8');
	im(2, fbp, 'FBP', climp), cbar
prompt
end

if ~isvar('fw.fbp'), printm 'fw.fbp'
	ox = [-9:9]; oy = ox;
	fw.fun = @(x, ii) fwhm2(x(ix(ii)+ox, iy(ii)+oy));
	for ii=1:numel(ix)
		fw.fbp(ii) = fw.fun(fbp, ii);
	end
	pr fw.fbp
end

if ~isvar('kappa'), printm 'kappa: try to make resolution approximately uniform'
	kappa = sqrt( div0(A' * wi, A' * ones(size(wi))) );
	im(4, kappa), cbar
prompt
end

% use local psf to help select beta
if ~isvar('R1'), printm 'R1, R2'
	f.l2b = 0;
	W = Gdiag(wi);
	R1 = Reg1(ig.mask * kappa(end/2+1,end/2+1), 'beta', 2^f.l2b); % usual
%	qpwls_psf(A, R1, 1, ig.mask, W, 'loop', 1); % choose beta
	R2 = Reg1(kappa, 'beta', 2^f.l2b); % kappa
	qpwls_psf(A, R2, 1, ig.mask, W, 'loop', 1); % choose beta
prompt
end

%	psf = qpwls_psf(A, R1, 1, ig.mask, W);
%	init = conv2(psf, xp, 'same');
	init = fbp;
	im(init, climp)

if ~isvar('xpwls1'), printm 'pwls1'
	f.niter = 300;
	xpwls1 = pwls_pcg1(init(ig.mask), A, W, yp(:), R1, 'niter', f.niter);
	xpwls1 = ig.embed(xpwls1);
	im(3, xpwls1, 'PWLS R1'), cbar
end

if ~isvar('xpwls2'), printm 'pwls2'
	xpwls2 = pwls_pcg1(init(ig.mask), A, W, yp(:), R2, 'niter', f.niter);
	xpwls2 = ig.embed(xpwls2);
	im(4, xpwls2, 'PWLS R2'), cbar
end

if 1
	for ii=1:numel(ix)
		fw.pwls1(ii) = fw.fun(xpwls1, ii);
		fw.pwls2(ii) = fw.fun(xpwls2, ii);
	end
	pr fw.fbp
	pr fw.pwls1
	pr fw.pwls2
end

if 1
	climp = [-0.004 0.05];
%	climp = []
	im plc 1 3
	ax = [18 110 20 120];
	im(1, fbp, 'FBP', climp), axis(ax)%, cbar h
	fun = @(fw, ii) text(ix(ii), iy(ii)+15, ...
		sprintf('%3.1f', fw(ii)), 'horiz', 'center', 'color', 'green');
	for ii=1:numel(ix), fun(fw.fbp, ii); end
	im(2, xpwls1, 'PWLS standard R', climp), axis(ax)%, cbar h
	for ii=1:numel(ix), fun(fw.pwls1, ii); end
	im(3, xpwls2, 'PWLS modified R', climp), axis(ax)%, cbar h
	for ii=1:numel(ix), fun(fw.pwls2, ii); end

%	ir_savefig eps_c ir_example_ct_lir1_fwhm
end
