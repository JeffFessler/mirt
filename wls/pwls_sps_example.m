%% pwls_sqs_example.m
%
% Complete example m-file illustrating algorithms for PWLS cost function
% (nonquadratically penalized weighted least squares).
% Algorithms: SQS (separable quadratic surrogates),
% SQS-OS (order subsets verions) and CG (conjugate gradents)
%
% Copyright 2001-8-16, Jeff Fessler, The University of Michigan

%% generate data
if ~isvar('yi'), printm 'setup pwls_sqs_example'
	f.count = 2e5;
	em_wls_test_setup
	W = diag_sp(wi(:));
prompt
end


%% regularizer
if ~isvar('Rn'), printm 'make R'
	f.l2b_n = 12; % log_2(beta)
	f.l2b_q = 9;

	% kappa terms from fessler:96:srp paper
	kappa = sqrt( (G' * wi(:)) ./ (G' * ones(size(wi(:)))) );
	kappa = ig.embed(kappa);
	im(8, kappa, 'kappa'), cbar

	Rq = Robject(kappa, 'potential', 'quad', 'beta', 2^f.l2b_q);

	% Play with regularization parameter (beta) to look at PSF.
	% I manually adjusted f.l2b above until I got fwhm=1.2 pixels or so.
	if 1
		psf = qpwls_psf(G, Rq, 1, ig.mask, W);
		printm('fwhm = %g', fwhm2(psf))
		im(7, psf, 'PSF'), cbar
		clear psf
	prompt
	end

	f.delta = 0.2;		% this depends on your units!
	f.pot = 'huber';	% nonquadratic, edge-preserving penalty function
	Rn = Robject(kappa, 'type_denom', 'matlab', ...
		'potential', f.pot, 'beta', 2^f.l2b_n, 'delta', f.delta);
end


%% Unconstrained conjugate gradient with quadratic penalty
if ~isvar('xpcg'), printm 'PCG'
	f.niter = 16;
	%xinit = ones(size(ig.mask));	% uniform initial image
	xinit = max(xfbp,0);		% FBP initial image

	xpcg = qpwls_pcg(xinit(ig.mask), G, W, yi(:), 0, ...
		Rq.C, 'circ0', f.niter, ig.mask);
	xpcg = ig.embed(xpcg);
	im clf, im(xpcg, 'QPWLS-CG'), cbar
prompt
end


%% SPS=SQS iterations
if ~isvar('xsqs'), printm 'do SQS'
	f.sqs_niter = 25;
	xsqs = pwls_sps_os(xinit(ig.mask), yi, wi, G, Rn, ...
		f.sqs_niter, inf, [], [], 1, 0);
	xsqs = ig.embed(xsqs);
	im clf, im(xsqs, 'PWLS-SQS'), cbar
prompt
end


%% SQS-OS iterations
% fix: replace with incremental version
if ~isvar('Gb'), printm 'make Gblock'
	f.nsubset = 8;
	Gb = Gblock(G, f.nsubset);
	clear tmp
end
if ~isvar('xsos'), printm 'SQS-OS'
	xsos = pwls_sps_os(xinit(ig.mask), yi, wi, Gb, Rn, ...
		f.sqs_niter, inf, [], [], 1, 0);
	xsos = ig.embed(xsos);
	im clf, im(xsos, 'PWLS-SQS-OS'), cbar
prompt
end


%% look at final results
if 1
	xpcgn = max(xpcg(:,:,end),0); % final nonnegativity for CG
	xfbpn = max(xfbp,0);
	im plc 2 5
	clim = [0 9];
	elim = [-3 3];
	im(1, xtrue, 'xtrue', clim), cbar horiz
	im(2, xfbpn, 'FBP', clim), cbar horiz
	im(3, xpcgn, 'QPWLS-PCG (+)', clim), cbar horiz
	im(4, xsqs(:,:,end), 'PWLS-SQS', clim), cbar horiz
	im(5, xsos(:,:,end), 'PWLS-SQS-OS', clim), cbar horiz
	efbp = xfbpn - xtrue;
	epcg = xpcgn - xtrue;
	esqs = xsqs(:,:,end) - xtrue;
	esos = xsos(:,:,end) - xtrue;
	nfbp = norm(efbp(:)) / norm(xtrue(:)) * 100;
	npcg = norm(epcg(:)) / norm(xtrue(:)) * 100;
	nsqs = norm(esqs(:)) / norm(xtrue(:)) * 100;
	nsos = norm(esos(:)) / norm(xtrue(:)) * 100;
	im(7, efbp, elim), cbar horiz, titlef('NRMSE %4.1f\%%', nfbp)
	im(8, epcg, elim), cbar horiz, titlef('NRMSE %4.1f\%%', npcg)
	im(9, esqs, elim), cbar horiz, titlef('NRMSE %4.1f\%%', nsqs)
	im(10, esos, elim), cbar horiz, titlef('NRMSE %4.1f\%%', nsos)
end
