% pwls_example.m
% example of how to use PWLS algorithms for PET reconstruction, such as PWLS-PCG
% Copyright 2003-11-30, Jeff Fessler, University of Michigan

% generate data
if ~isvar('yi'), printm 'setup'
	em_wls_test_setup
%	wi = ones(size(wi)); warning 'uniform wi' % to test circulant precon
	W = diag_sp(wi(:));
prompt
end

% regularization
if ~isvar('R'), printm 'R'
	f.l2b = 9;
	f.delta = 1;

	if 1
		Rq = Reg1(ig.mask, 'beta', 2^f.l2b);
		psf = qpwls_psf(G, Rq.C, 1, ig.mask);
		im(7, psf, 'PSF'), cbar
	end, clear Rq

	kappa = sqrt( div0(G' * wi, G' * sg.ones) );
	im(8, kappa, 'kappa'), cbar

	R = Reg1(kappa, 'type_denom', 'matlab', ...
		'pot_arg', {'hyper3', f.delta}, 'beta', 2^f.l2b);

	clear xos xiot xcg
prompt
end

f.niter = 20;

if ~isvar('xinit')
	%xinit = ones(size(xtrue));		% uniform initial image
	%xinit = xfbp;				% FBP initial image
	if ir_has_imfilter
		xinit = imfilter(xfbp, double(psf)); % smoothed FBP
	else
		xinit = convn(xfbp, double(psf), 'same'); % smoothed FBP
	end
	im(9, xinit, 'init'), cbar
prompt
end


if ~isvar('Ab'), printm 'do Ab'
	f.nsubset = 40;
	Ab = Gblock(G, f.nsubset);
prompt
end


% OS-SPS iterations (unconstrained)
if ~isvar('xos'), printm 'do os'
	xlims = [-inf inf]; % unconstrained
	[xos tim.os] = pwls_sqs_os(xinit(ig.mask), Ab, yi, R, 'wi', wi, ...
			'niter', f.niter, 'pixmax', xlims, 'isave', 'all', 'chat', 1);
	xos = ig.embed(xos);
	im clf, im(xos, 'xos'), cbar
prompt
end


% PWLS-IOT iterations (unconstrained)
if 0 || ~isvar('xiot'), printm 'do iot'
	[xiot tim.iot] = pl_iot(xinit(ig.mask), Ab, {yi, wi}, R, ...
			'dercurv', 'wls', ...
			'riter', 1, ...
			'os', 5, ... % f.niter, ...
			'niter', f.niter, 'isave', 'all', ...
			'pixmin', xlims(1), 'pixmax', xlims(2), ...
			'chat', 0);
	xiot = ig.embed(xiot);
	im clf, im(xiot, 'xiot'), cbar
	minmax(xiot-xos)
prompt
end


% CG iterations
if ~isvar('xcg'), printm 'xcg'
	[xcg tim.cg] = pwls_pcg1(xinit(ig.mask), G, W, yi(:), R, ...
			'niter', f.niter, 'isave', 'all');
	xcg = ig.embed(xcg);
	im clf, im(xcg, 'CG'), cbar
prompt
end

% compare cost
if 1, printm 'cost plots'
	cost.fun = @(x) pwls_cost(x, G, W, yi(:), R, ig.mask);
	cost.cg	 = cost.fun(xcg);
	cost.os	 = cost.fun(xos);
	cost.iot = cost.fun(xiot);
	ii = 0:f.niter;
	if im
		clf, subplot(211)
		plot(ii, cost.os, 'y-o', ii, cost.cg, 'g-x', ii, cost.iot, 'c-+')
		xlabel 'iteration', ylabel 'cost'
		legend('OS-SPS', 'CG', 'IOT')

		subplot(212)
		plot([0; tim.os], cost.os, 'y-o', [0; tim.cg(:,3)], cost.cg, 'g-x', ...
			tim.iot, cost.iot, 'c-+')
		xlabel 'time', ylabel 'cost'
		legend('OS-SPS', 'CG', 'IOT')
	prompt
	end
%	minmax(diff(tim.os))
%	minmax(diff(tim.cg))
end

	pr 'nrms(col(xos(:,:,end)), xtrue(:))'
	pr 'nrms(col(xcg(:,:,end)), xtrue(:))'
	pr 'nrms(col(xiot(:,:,end)), xtrue(:))'

if 1, printm 'images'
	clim = [0 6];
	im clf, im([xtrue, xfbp; xos(:,:,end), xiot(:,:,end); ...
		xcg(:,:,end), xinit], clim), cbar
	im plc 2 3
	im(1, xtrue, clim)
	im(2, xfbp, clim, 'FBP')
	im(3, xinit(:,:,end), clim, 'Init')
	im(4, xos(:,:,end), clim, 'OS')
	im(5, xiot(:,:,end), clim, 'IOT')
	im(6, xcg(:,:,end), clim, 'CG')
end
