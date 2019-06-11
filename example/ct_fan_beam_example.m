% ct_fan_beam_example.m
% compare FBP and iterative reconstruction for a 2D fan-beam CT problem
% Copyright 2009-10-05, Jeff Fessler, University of Michigan

if ~isvar('A'), printm 'setup geometry, image, sinogram'
	f.down = 4;
	f.do_sparse = 1; % set to one for sparse angle case
	ig = image_geom('nx', 512, 'fov', 50, 'down', f.down);
	ig.mask = ig.circ > 0;
	sg = sino_geom('ge1', 'units', 'cm', 'strip_width', 'd', ...
		'down', f.down);
	if f.do_sparse, sg.na = 50; end % # views in sparse angle case

	% system object
	if has_mex_jf
		A = Gtomo2_dscmex(sg, ig);
	else
		A = Gtomo_nufft_new(sg, ig);
	end

	% read image
	ddir = path_find_dir([filesep 'data']);
	xtrue256 = fld_read([ddir filesep 'ncat,256,slice,140,ct,x100.fld']);
	xtrue256 = single(xtrue256) / 200 * 0.4; % convert to 1/cm units 

	if 1 % more realistic sinogram from finer image, avoid "inverse crime"
		ig_big = image_geom('nx', 512, 'fov', ig.fov, 'down', 2);
		if has_mex_jf
			Abig = Gtomo2_dscmex(sg, ig_big);
		else
			Abig = Gtomo_nufft_new(sg, ig_big);
		end
		sino_true = Abig * xtrue256;
	end
	xtrue = downsample2(xtrue256, 2);

	im plc 2 2
	clim = [0 0.4];
	im(1, xtrue, 'x', clim), cbar
	xlabelf('units: 1 / %s', sg.units)
	im(2, sino_true, 'sino'), cbar

	clear ddir ig_big Abig
prompt
end


if ~isvar('sino'), printm 'noisy fan-beam data'
	I0 = 1e5; % incident photons; decrease this for "low dose" scans
	rng(0)
	% transmission data:
	yi = poisson(I0 * exp(-sino_true), 0, 'factor', 0.4); % poisson noise
	if any(yi(:) == 0)
		warn('%d of %d values are 0 in sinogram!', ...
			sum(yi(:)==0), length(yi(:)));
	end
	sino = log(I0 ./ max(yi,1)); % noisy fan-beam sinogram
	im(4, sino, 'noisy sino'), cbar
prompt
end

xl = @(x) xlabelf('RMSE = %.3f / %s', rms(col(x - xtrue)), sg.units);

if ~isvar('fbp'), printm 'fbp 2d fan-beam reconstruction'
	tmp = fbp2(sg, ig);
	fbp = fbp2(sino, tmp);
	im(3, fbp, 'FBP', clim), cbar
	xl(fbp)
prompt
end

if f.do_sparse && ~isvar('fbpw') % sparse angle
	tmp = fbp2(sg, ig);
	fbpw = fbp2(sino, tmp, 'window', 'hanning,0.6');
	im(4, fbpw, 'FBP Hanning', clim), cbar
	xl(fbpw)
prompt
end


if ~isvar('kappa'), printm 'kappa: try to make resolution approximately uniform'
	wi = yi; % will give 0 weight to any ray where yi=0!
	kappa = sqrt( div0(A' * wi, A' * ones(size(wi))) );
	im(4, kappa), cbar
prompt
end


% use local psf to help select beta
if ~isvar('R'), printm 'R'
	f.l2b = 8; % maybe a bit too big, but ok for now
	f.delta = 0.01;
	if f.do_sparse
		f.l2b = 4; % sparse angle
	end
%	f.pot_arg = {'lange3', f.delta}; % todo: why not as sharp as hyper3?
	f.pot_arg = {'hyper3', f.delta};
	R = Reg1(kappa, 'beta', 2^f.l2b, 'pot_arg', f.pot_arg);
%	qpwls_psf(A, R, 1, ig.mask, Gdiag(wi), 'loop', 1); % use this to choose beta
end


if ~isvar('xpwls'), printm 'iterative reconstruction'
	if has_mex_jf
		f.niter = 20;
		Ab = Gblock(A, 41); % 41 subsets
		xpwls = pwls_sqs_os(fbp(ig.mask), Ab, sino, R, 'wi', wi, 'niter', f.niter);
	else
		f.niter = 90;
		xpwls = pwls_pcg1(fbp(ig.mask), A, Gdiag(wi), sino(:), R, 'niter', f.niter);
	end
	xpwls = ig.embed(xpwls);
	im(4, xpwls, 'PWLS'), cbar
	xl(xpwls)
%prompt
end


if f.do_sparse % sparse angle figure for book chapter
	im plc 2 2
	im(1, xtrue, 'x', clim), cbar
	xlabelf('units: 1 / %s', sg.units)
	im(2, xpwls, 'PWLS', clim), cbar
	xl(xpwls)
	im(3, fbp, 'FBP Ramp', clim), cbar
	xl(fbp)
	im(4, fbpw, 'FBP Hanning', clim), cbar
	xl(fbpw)
%	ir_savefig ct_fan_beam_example1
prompt
end
