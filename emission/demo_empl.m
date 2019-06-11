% demo_empl.m
% interactive demo for penalized-likelihood emission image reconstruction
%
% Copyright 2010-03-22, Jeff Fessler, University of Michigan


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


% construct gui components
if 1 || ~isvar('h_l2b')
	im clf
	im(xinit, f.clim, 'Init')
	h_l2b = uicontrol('style', 'slider', 'string', 'value', ...
		'units', 'normalized', 'position', [0.1 0.0 0.8 0.03], ...
		'value', 0);

	h_del = uicontrol('style', 'slider', 'string', 'value', ...
		'units', 'normalized', 'position', [0.1 0.0 0.8 0.06], ...
		'value', 1);
prompt
end


% run regularized case for comparison
if 1 || ~isvar('xh')

	save_l2b = nan;
	save_del = nan;

	f.fun_l2b = @(x) -4 + 6 * x;
	f.fun_del = @(x) 0.01 + 1.5 * x;

	xh = xinit;
	aj = Gb' * ci(:) / f.nblock;

	while (1)
		l2b = f.fun_l2b( get(h_l2b, 'value') );
		del = f.fun_del( get(h_del, 'value') );
		if (l2b ~= save_l2b || del ~= save_del)
			iter = 0;
			save_l2b = l2b;
			save_del = del;
	%		printm 'build R'
			R = Reg1(ig.mask, 'type_denom', 'matlab', ...
				'pot_arg', {'hyper3', del}, 'beta', 2^l2b);
		end

		f.niter_pl = 2;
		xh = epl_os_emdp(xh(ig.mask), Gb, yi, ci, ri, R, ...
			'niter', f.niter_pl, 'aj', aj, ...
			'userfun', [], 'isave', 'last');
		xh = ig.embed(xh);

		iter = iter + f.niter_pl - 1;
		im(xh, f.clim)
		titlef('iter=%4d log_2(beta) = %g', iter, l2b)
		texts(0.05, 0.04, sprintf('\\delta = %g', del))
		texts(0.05, 0.08, sprintf('log_2(\\beta) = %g', l2b))
		drawnow
	end
end

return
% old below here ......................

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
