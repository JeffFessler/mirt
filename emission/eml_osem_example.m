% eml_osem_example.m
%
% a complete example m-file illustrating E-ML-OSEM.
% showing faster "convergence" compared to E-ML-EM
%
% Copyright 1998-1, Jeff Fessler, University of Michigan


% generate data
if ~isvar('yi'), printm 'set up'
	f.count = 1e6;
	im clf, em_test_setup
	f.clim = [0 6];
prompt
end


% ML-EM iterations
if ~isvar('xmlem'), printm 'E-ML-EM'
	f.niter = 100;
	xinit = ig.ones; % uniform initial image

	im clf, im(xinit, f.clim, 'init')
	fun = @jf_show_iter;
	fun('mask', ig.mask, 'clim', f.clim, 'pause', 0); % initialize

	xmlem = eml_em(xinit(ig.mask), G, yi(:), ci(:), ri(:), ...
		'niter', f.niter, 'isave', 'all', 'userfun', fun);
	xmlem = ig.embed(xmlem);

	im plc 1 2
	im(1, xtrue, f.clim), cbar
	im(2, xmlem(:,:,end), 'ML-EM', f.clim), cbar
	printm('Done running ML-EM')
prompt
end


% block object for block iterations
if ~isvar('Gb'), printm 'Gb'
	f.nblock = 5;
	Gb = Gblock(G, f.nblock);
end


% OS-EM iterations
if ~isvar('xosem'), printm 'E-ML-OS-EM'
	im clf, im(xinit, f.clim, 'init')
	xosem = eml_osem(xinit(ig.mask), Gb, yi, ci, ri, ...
		'niter', f.niter, 'isave', 'all', 'precon', 'classic', ...
		'userfun', fun);
	xosem = ig.embed(xosem);

	im plc 1 3
	im(1, xtrue, f.clim), cbar
	im(2, xmlem(:,:,end), f.clim, 'ML-EM'), cbar
	im(3, xosem(:,:,end), f.clim, 'E-ML-OS-EM classic'), cbar
prompt
end


if 0, printm 'E-ML-OS-EM fast'
	xfast = eml_osem(xinit(ig.mask), Gb, yi, ci, ri, ...
		'niter', f.niter-1, 'precon', 'fast');
	xfast = ig.embed(xfast);

	im(3, xfast, 'E-ML-OS-EM fast'), cbar
	im(4, xfast-xosem, 'fast-classic'), cbar
prompt
end


% E-ML-INC-EM-3 iterations
if ~isvar('xiem3'), printm 'E-ML-INC-EM'

	xiem1 = eml_inc_em(xinit(ig.mask), Gb, yi, ci, ri, ...
		'niter', f.niter, 'hds', 1, 'isave', 'all');
	xiem1 = ig.embed(xiem1);

	xiem3 = eml_inc_em(xinit(ig.mask), Gb, yi, ci, ri, ...
		'niter', f.niter, 'hds', 3, 'isave', 'all');
	xiem3 = ig.embed(xiem3);

	im plc 1 2
	im(1, xiem1(:,:,end), f.clim, 'E-ML-INC-EM-1'), cbar
	im(2, xiem3(:,:,end), f.clim, 'E-ML-INC-EM-3'), cbar
	printm 'Done running E-ML-INC-EM'
prompt
end


if ~isvar('iplot') % plot likelihood to show acceleration
	iplot = 0:10;
	lfun = @(x) eql_obj(x(:,:,iplot+1), G, yi(:), ci(:), ri(:), [], ig.mask);
	like.mlem = lfun(xmlem);
	like.osem = lfun(xosem);
	like.iem1 = lfun(xiem1);
	like.iem3 = lfun(xiem3);
	if im
		im clf, plot(	...
			iplot, like.osem, 'c-+', ...
			iplot, like.iem1, 'g*-', ...
			iplot, like.iem3, 'm*-', ...
			iplot, like.mlem, 'yx-', ...
			f.nblock * iplot, like.osem, 'ro')
		axisx(0, max(iplot))
		legend('OS-EM', 'INC-EM-1', 'INC-EM-3', 'ML-EM', ...
			sprintf('OSEM/%d', f.nblock), 4)
		title 'ML-EM vs ML-OSEM convergence rate'
		xlabel iteration, ylabel likelihood
	end
prompt
end


% run regularized case for comparison
if ~isvar('R'), printm 'build R'
	f.l2b = -2;
	f.delta = 0.5;
	R = Reg1(ig.mask, 'type_denom', 'matlab', ...
		'pot_arg', {'hyper3', f.delta}, 'beta', 2^f.l2b);
%prompt
end


if ~isvar('xh'), printm 'E-PL-OS-EMDP'
	f.niter_pl = 31;
	im clf
	xh = epl_os_emdp(xinit(ig.mask), Gb, yi, ci, ri, R, ...
		'niter', f.niter_pl, 'userfun', fun, 'isave', 'last');
	xh = ig.embed(xh);
	im plc 1 3
	im(1, xtrue, f.clim), cbar
	im(2, xmlem(:,:,end), f.clim), cbar
	im(3, xh, 'PL-OS-EMDP', f.clim), cbar
%prompt
end

if ~isvar('elim') % compare pics of ml-osem vs pl
	elim = [-1 1]*2.5;
	im plc 2 3
	im(1, xtrue, f.clim, 'true'), cbar
	im(2, xosem(:,:,end), f.clim, 'ML-OSEM'), cbar
	im(3, xh(:,:,end), f.clim, 'PL-OS-EMDP'), cbar
	im(5, xtrue-xosem(:,:,end), elim), cbar
	titlef('nrms err %g%%', 100*nrms(col(xosem(:,:,end)), xtrue(:)))
	im(6, xtrue-xh(:,:,end), elim), cbar
	titlef('nrms err %g%%', 100*nrms(col(xh(:,:,end)), xtrue(:)))
prompt
end

if 1
	im clf
	im([xtrue; xmlem(:,:,end); xh(:,:,end)], f.clim, 'True | ML | PL')
prompt
end

if ~isvar('xcake') % illustrate space-variant PSF
	[px py] = ndgrid(3:7:ig.nx, 3:7:ig.ny);
	xcake = ig.zeros;
	xcake(px(:), py(:)) = 1;
	xcake = xcake .* ig.mask;
	im(xcake)

	if ~isvar('xh2'), printm 'E-PL-OS-EMDP'
		im clf
		frac = 0.4;
		ytmp = yi + frac * ci .* (G * xcake);
		xh2 = epl_os_emdp(xinit(ig.mask), Gb, ytmp, ci, ri, R, ...
			'niter', f.niter_pl, 'userfun', fun, 'isave', 'last');
		xh2 = ig.embed(xh2);
	end
%prompt
end
	im plc 2 2
	im(1, xh, 'PL Original', f.clim), cbar
	im(2, xcake), cbar
	im(3, xh2, 'PL with point perturbations', f.clim), cbar
	im(4, (xh2 - xh) / frac, 'PSF', [-0.1 0.5]), cbar
