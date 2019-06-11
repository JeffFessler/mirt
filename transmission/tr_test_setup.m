% tr_test_setup.m
%
% create sample image, system matrix, and sinograms for examples
% and testing of Poisson transmission maximum likelihood (ML) algorithms
% creates: ig sg xtrue G bi ri yi
%
% Copyright Apr 2000, Jeff Fessler, University of Michigan

%
% true transmission image
% use data/do,attn,zubal script to 'generate' data if necessary
%
if ~isvar('xtrue'), printm 'xtrue'
	ig = image_geom('nx', 64, 'ny', 60, 'fov', 500);

	xtrue = read_zubal_attn('nx', ig.nx, 'ny', ig.ny);

	im pl 3 3
	im(1, xtrue, 'true attenuation image'), cbar

	% reconstruction mask (which pixels do we estimate?)
	ig.mask = conv2(double6(xtrue > 0), ones(3), 'same') > 0;
	im(4, ig.mask, 'support mask')
end

%
% system matrix G, scaled by "dx" so it has "line length" units
%
if ~isvar('G'), printm 'G'
	sg = sino_geom('par', 'nb', ig.nx+2, 'dr', 528 / (ig.nx+2));

	% simple strip-integral system model
	G = Gtomo2_strip(sg, ig, 'single', 1);

	if isvar('f.wtf') && has_mex_jf
		if exist(f.wtf), delete(f.wtf), end
		wtf_write(f.wtf, G, ig.nx, ig.ny, sg.nb, sg.na);
	end
	if isvar('f.wtr') && has_mex_jf && has_aspire
		if exist(f.wtr, 'file'), delete(f.wtr), end
		os_run(sprintf('wt -chat 0 col2row %s %s', f.wtr, f.wtf))
	end
end


%
% noisy measurements
%
if ~isvar('yi')
	proj = G * xtrue;
	printm('Maximum line integral = %g', max(proj(:)))
	if ~isvar('f.count'), f.count = 1e6; end
	bi = f.count / sum(exp(-proj(:))) * sg.ones;
	bi = dsingle(bi);
	ytrue = bi .* exp(-proj);
	if ~isvar('f.randpercent'), f.randpercent = 10; end
	ri = f.randpercent / 100 * mean(ytrue(:)) * sg.ones;
	ri = dsingle(ri);
	rng(0)
	yi = poisson(ytrue + ri);

	im(2, ytrue, 'ytrue: true projections'), cbar
	im(3, yi, 'yi: noisy projections'), cbar
	clear ytrue proj
end

%
% FBP reconstruction
%
if ~isvar('xfbp')
	xfbp = tr_fbp(sg, ig, yi, bi, ri);
	xfbp = max(xfbp, 0);
	im(5, xfbp, 'FBP Reconstruction'), cbar
end

if isvar('f.yi')
	fld_write(f.yi, yi, 'check', 0)
end
if isvar('f.bi')
	fld_write(f.bi, bi, 'check', 0)
end
if isvar('f.ri')
	fld_write(f.ri, ri, 'check', 0)
end
if isvar('f.mask')
	fld_write(f.mask, ig.mask, 'check', 0)
end
