% qpwls_qn_test.m
% Compare Quasi-Newton (QN) and SPS algorithms for QPWLS problem.
% Also shows the benefit of initializing QN with good diagonal Hessian approx.
% Copyright July 2000, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'), printm 'setup qpwls_qn_test'
	em_wls_test_setup
	W = diag_sp(wi(:));
prompt
end

if ~isvar('R'), printm 'make R'
	f.l2b = -2;
	R = Robject(ig.mask, 'type_denom', 'matlab', 'beta', 2^f.l2b);
prompt
end


if ~isvar('xinit')
	f.niter = 12;
	xinit = xfbp;	% warning: no nonnegativity - QN unconstrained!
end

if ~isvar('P0')
	P0 = G' * (W * sum(G')') + R.denom(R, 0);
	P0 = diag_sp(1 ./ P0(:));
end

%
% PCG
%
if ~isvar('xpcg'), printm 'do pcg'
	xpcg = qpwls_pcg(xinit(ig.mask), G, W, yi, 0, R.C, P0, 1+f.niter);
	xpcg = ig.embed(xpcg);
	im clf, im(xpcg, 'QPWLS-PCG iterations')
prompt
end


%
% test new PCG
%
if ~isvar('xnew')
	if ~isvar('xnew'), printm 'do pcg new'
		xnew = pl_pcg_qs_ls(xinit(ig.mask), G, {yi(:), wi(:)}, ...
			@wls_dercurv, R, 'precon', P0, 'niter', f.niter, ...
			'isave', 'all');
		xnew = ig.embed(xnew);
	end
	max_percent_diff(xpcg, xnew)
prompt
end


%
% QN
%
if ~isvar('xqnp'), printm 'do qnp'
	xqnp = qpwls_qn(xinit(ig.mask), G, W, yi, R.C, P0, 1+f.niter);
	xqnp = ig.embed(xqnp);
	im clf, im(xqnp, 'QPWLS-QNP iterations')
prompt
end
if ~isvar('xqnu'), printm 'do qnu'
	xqnu = qpwls_qn(xinit(ig.mask), G, W, yi, R.C, 1, 1+f.niter);
	xqnu = ig.embed(xqnu);
	im clf, im(xqnu, 'QPWLS-QNU iterations')
prompt
end

%
% SPS
%
if ~isvar('xsps'), printm 'do sps'
	xsps = pwls_sps_os(xinit(ig.mask), yi, wi, G, R, ...
		1+f.niter, [-inf inf], [], [], 1, 0); % disable nonnegativity
	xsps = ig.embed(xsps);
	im clf, im(xsps, 'QPWLS-SPS iterations')
prompt
end


if 1
	im pl 2 2, clim = [0 8];
	im(1, xsps, 'SPS matlab', clim), cbar
	im(2, xqnp, 'QNP matlab', clim), cbar
	im(3, xpcg, 'PCG matlab', clim), cbar
	im(4, xqnu, 'QNU matlab', clim), cbar
%	im(223, (xart2-xmat)/max(xmat(:)), '(aspire-matlab)/max(matlab)')

	cost = @(x) pwls_cost(x, G, W, yi(:), R, ig.mask);
	t1 = cost(xsps);
	t2 = cost(xpcg);
	t3 = cost(xqnp);
	t4 = cost(xqnu);
	if im
		subplot(224)
		plot(	0:f.niter, t1-t1(1), '-o', ...
			0:f.niter, t2-t1(1), '-+', ...
			0:f.niter, t3-t1(1), '-x', ...
			0:f.niter, t4-t1(1), '-v')
		xlabel iteration, ylabel '\Phi change',
		ir_legend({'SPS', 'PCG', 'QNP', 'QNU'})
		title('QPWLS: QN vs SPS vs PCG')
%		axis([0 min(f.niter, 20) -18000 -15000])
		axisy(-5000, 500)
	end

%	ir_savefig 'fig_qpwls_qn_vs_sps'
end
