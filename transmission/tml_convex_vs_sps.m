% tml_convex_vs_sps.m
% Compare T-ML-Convex-{PS,Log} and T-ML-SPS,
% all of which are monotone algorithms!
% Copyright 2001-10-16, Jeff Fessler, University of Michigan


% generate data
if ~isvar('yi')
	f.randpercent = 10;
	tr_test_setup
	xinit = 0.01+0*xfbp; % uniform initial image
%	xinit = xfbp;
end

f.niter = 15;
f.pixmax = 0.012;


% sps
if ~isvar('xsps'), printm 'xsps'
	cpu etic
	xsps = tpl_os_sps(xinit(ig.mask), G, yi, bi, ri, [], ...
		f.niter+1, f.pixmax, 'oc');
	cpu etoc 'sps time'
	xsps = ig.embed(xsps);
	im clf, im(221, xsps, 'T-ML-SPS')
	lsps = tpl_obj(xsps, G, yi(:), bi(:), ri(:), [], ig.mask);
prompt
end


% convex-ps
if ~isvar('xc.ps'), printm 'xc.ps'
	cpu etic
	xc.ps = tml_convex(xinit(ig.mask), G, yi, bi, ri, f.niter+1, ...
		f.pixmax, 'ps');
	cpu etoc 'convex time'
	xc.ps = ig.embed(xc.ps);
	im(222, xc.ps, 'T-ML-Convex-PS')
	lc.ps = tpl_obj(xc.ps, G, yi(:), bi(:), ri(:), [], ig.mask);
prompt
end


% convex-log1
if ~isvar('xc.l1'), printm 'xc.l1'
	cpu etic
	xc.l1 = tml_convex(xinit(ig.mask), G, yi, bi, ri, f.niter+1, ...
		f.pixmax, 'log1');
	cpu etoc 'convex log1 time'
	xc.l1 = ig.embed(xc.l1);
	im(223, xc.l1, 'T-ML-Convex-Log1')
	lc.l1 = tpl_obj(xc.l1, G, yi(:), bi(:), ri(:), [], ig.mask);
prompt
end


% convex-log2
if ~isvar('xc.l2'), printm 'xc.l2'
	cpu etic
	xc.l2 = tml_convex(xinit(ig.mask), G, yi, bi, ri, f.niter+1, f.pixmax, 'log2');
	cpu etoc 'convex log2 time'
	xc.l2 = ig.embed(xc.l2);
	im(223, xc.l2, 'T-ML-Convex-Log2')
	lc.l2 = tpl_obj(xc.l2, G, yi(:), bi(:), ri(:), [], ig.mask);
prompt
end


% figure comparing convergence
if im
	im clf
	plot(	...
		0:f.niter, lsps-lsps(1), '.-', ...
		0:f.niter, lc.ps-lsps(1), '-o', ...
		0:f.niter, lc.l1-lsps(1), '-^', ...
		0:f.niter, lc.l2-lsps(1), '-d')
%	axisy([0.6 1.9] * 1e5)
	ir_legend({'SPS', 'Convex-PS', 'Convex-Log1', 'Convex-Log2'})
	xlabel 'Iteration', ylabel 'Likelihood Increase'
	title 'Transmission ML via Monotonic Algorithms'
	axes('position', [0.6 0.4 0.2 0.2])
	im('notick', xtrue), title 'Object'
	colormap(1-gray(256)), cbar
%	ir_savefig 'tml_convex_vs_sps'
end


% now look at nonmonotone algorithms
if 0
	% SPS-nr
	xspsnr = tpl_os_sps(xinit(ig.mask), G, yi, bi, ri, [], ...
		f.niter+1, f.pixmax, 'nr');
	xspsnr = ig.embed(xspsnr);
	lspsnr = tpl_obj(xspsnr, G, yi(:), bi(:), ri(:), [], ig.mask);

	% convex-nr2
	xc.nr2 = tml_convex(xinit(ig.mask), G, yi, bi, ri, f.niter+1, f.pixmax, 'nr2');
	xc.nr2 = ig.embed(xc.nr2);
	lc.nr2 = tpl_obj(xc.nr2, G, yi(:), bi(:), ri(:), [], ig.mask);

	im clf
	plot(	...
		0:f.niter, lspsnr-lsps(1), '-^', ...
		0:f.niter, lc.nr2-lsps(1), '-x', ...
		0:f.niter, lsps-lsps(1), '-s')
	ir_legend({'SPS-NR', 'Convex-NR', 'SPS'})
	xlabel Iteration, ylabel 'Likelihood Increase'
	title 'Transmission ML via Non-Monotonic Algorithms'
end
