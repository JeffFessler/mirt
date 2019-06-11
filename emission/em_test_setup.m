% em_test_setup.m
%
% create sample image, system matrix, and sinograms for examples
% and testing of Poisson emission maximum likelihood (ML) algorithms
% creates: ig sg xtrue G proj ci ytrue ri yi
%
% Copyright Jan 1998, Jeff Fessler, University of Michigan

% true emission image
if ~isvar('xtrue'), printm 'xtrue'
	if ~isvar('ig')
		ig = image_geom('nx', 64, 'ny', 60, 'fov', 500);
	end
	xtrue = read_zubal_emis('nx', ig.nx, 'ny', ig.ny);
	mumap = read_zubal_attn('nx', ig.nx, 'ny', ig.ny);
	im plc 2 3
	im(1, xtrue, 'emission image'), cbar
	im(2, mumap, 'attenuation map'), cbar

	% reconstruction mask (which pixels do we estimate?)
	ig.mask = ig.circ(220, 180) > 0;
	im(3, ig.mask + xtrue, 'support mask + xtrue')
end


% system matrix G
if ~isvar('G'), printm 'system'
	sg = sino_geom('par', 'nb', ig.nx+2, 'na', ig.ny*3/2, ...
		'dr', 528 / (ig.nx+2));

	% simple strip-integral system model
	G = Gtomo2_strip(sg, ig, 'single', 1);

	if isvar('f.wtf') && has_mex_jf
		if exist(f.wtf), delete(f.wtf), end
		wtf_write(f.wtf, G, ig.nx, ig.ny, sg.nb, sg.na);
	end
	if isvar('f.wtr') && has_mex_jf && has_aspire
		if exist(f.wtr), delete(f.wtr), end
		os_run(sprintf('wt -chat 0 col2row %s %s', f.wtr, f.wtf))
	end
end


% noisy measurements
if ~isvar('yi'), printm 'data yi'
	proj = G * xtrue;
	li = G * mumap;
	printm('Maximum line integral = %g', max(li(:)))
	if ~isvar('f.count'), f.count = 1e5; end
	% detector efficiency variations per CTI 931 PET scanner
	ci = exp(0.3 * randn(size(proj)));
	ci = ci .* exp(-li);
	ci = f.count / sum(ci(:) .* proj(:)) * ci;
	ci = dsingle(ci);
	ytrue = ci .* proj;
	if ~isvar('f.randpercent')
		f.randpercent = 10;
	end
	ri = f.randpercent / 100 * mean(ytrue(:)) * sg.ones;
	ri = dsingle(ri);
	rng(0)
	yi = poisson(ytrue + ri);

	im(4, ytrue, 'ytrue: true projections'), cbar
	im(5, yi, 'yi: noisy projections'), cbar
	clear ytrue proj
end

% FBP reconstruction
if ~isvar('xfbp'), printm 'fbp'
	xfbp = em_fbp(sg, ig, yi, ci, ri);
	xfbp = max(xfbp, 0);
	im(6, xfbp, 'FBP Reconstruction'), cbar
end

% save to files if needed
if isvar('f.yi')
	fld_write(f.yi, yi, 'check', 0)
end
if isvar('f.ci')
	fld_write(f.ci, ci, 'check', 0)
end
if isvar('f.ri')
	fld_write(f.ri, ri, 'check', 0)
end
if isvar('f.mask')
	fld_write(f.mask, ig.mask, 'check', 0)
end
