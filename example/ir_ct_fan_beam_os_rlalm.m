% ir_ct_fan_beam_os_rlalm.m
% compare FBP and relaxed OS-LALM for a 2D fan-beam CT problem
% see ir_pwls_os_rlalm.m
% 2015-11-30 Hung Nien

if ~isvar('A'), printm 'setup geometry, image, sinogram'
	f.down = 4;
	ig = image_geom('nx', 512, 'fov', 50, 'down', f.down);
	ig.mask = ig.circ > 0;
	sg = sino_geom('ge1', 'units', 'cm', 'strip_width', 'd', ...
		'down', f.down);

	% system object
	if has_mex_jf
		A = Gtomo2_dscmex(sg, ig);
	else
		% Gblock(A) fails below for this so just abort here
		fail 'todo: requires mex file for forw/back_block()'
		A = Gtomo_nufft_new(sg, ig);
	end

	% read image
	ddir = path_find_dir([filesep 'data']);
	xtrue1024 = fld_read([ddir filesep 'ncat,1024,slice500,0to4.fld']);
	xtrue1024 = single(xtrue1024) / 10; % convert to 1/cm units 

	if 1 % more realistic sinogram from finer image, avoid "inverse crime"
		ig_big = image_geom('nx', 1024, 'fov', ig.fov);
		if has_mex_jf
			Abig = Gtomo2_dscmex(sg, ig_big);
		else
			Abig = Gtomo_nufft_new(sg, ig_big);
		end
		sino_true = Abig * xtrue1024;
	end
	xtrue = downsample2(xtrue1024, f.down*2);

	im plc 2 3
	clim = [0 0.4];
	im(1, xtrue, 'x', clim), cbar
	xlabelf('units: 1 / %s', sg.units)
	im(2, sino_true, 'sino'), cbar

	clear ddir ig_big Abig
drawnow
end


if ~isvar('sino'), printm 'noisy fan-beam data'
	rng(0)
	% transmission data:
	I0 = 1e4;
	yi = poisson(I0 * exp(-sino_true), 0, 'factor', 0.1);
	if any(yi(:) == 0)
		warn('%d of %d values are 0 in sinogram!', ...
			sum(yi(:)==0), length(yi(:)));
	end
	sino = log(I0 ./ max(yi,1)); % noisy fan-beam sinogram
	im(3, sino, 'noisy sino'), cbar
drawnow
end

xl = @(x) xlabelf('RMSE = %.3f / %s', rms(col(x - xtrue)), sg.units);

if ~isvar('fbp'), printm 'fbp 2d fan-beam reconstruction'
	tmp = fbp2(sg, ig);
	fbp = fbp2(sino, tmp, 'window', 'hanning,0.75');
	im(2, fbp, 'FBP', clim), cbar
	xl(fbp)
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
	f.l2b = 10; % maybe a bit too big, but ok for now
	f.delta = 0.001;
	f.pot_arg = {'hyper3', f.delta};
	R = Reg1(kappa, 'beta', 2^f.l2b, 'pot_arg', f.pot_arg);
end


% recon parameters
f.niter = 4;
f.nblock = 41; % 41 subsets
Ab = Gblock(A, f.nblock);

alpha = 1.999;

% rho = 1;	% ~ sqs
% rho = 0.05;	% 20x faster than sqs
rho = [];	% continuation

% OS-SQS
if ~isvar('xsqs'), printm 'iterative reconstruction - OS-SQS'
	xsqs = pwls_sqs_os(fbp(ig.mask), Ab, sino, R, 'wi', wi, ...
		'isave', 'all', 'niter', f.niter);
	xsqs = ig.embed(xsqs);
	im(3, xsqs(:,:,end), 'SQS', clim), cbar
	xl(xsqs(:,:,end))
prompt
end


% OS-LALM
if ~isvar('xlalm'), printm 'lalm recon'
	xlalm = ir_pwls_os_rlalm(fbp(ig.mask), Ab, sino, R, 'wi', wi, ...
		'alpha', 1, 'rho', rho, ...
		'isave', 'all', 'niter', f.niter);
	xlalm = ig.embed(xlalm);
	im(5, xlalm(:,:,end), 'LALM', clim), cbar
	xl(xlalm(:,:,end))
prompt
end


% relaxed OS-LALM
if ~isvar('xrla'), printm 'relaxed lalm recon'
	xrla = ir_pwls_os_rlalm(fbp(ig.mask), Ab, sino, R, 'wi', wi, ...
		'alpha', alpha, 'rho', rho, ...
		'isave', 'all', 'niter', f.niter);
	xrla = ig.embed(xrla);
	im(6, xrla(:,:,end), 'Relaxed LALM', clim), cbar
	xl(xrla(:,:,end))
prompt
end


if 1 && im
	if ~isvar('cost_sqs')
		cost_sqs = pwls_cost(xsqs, A, Gdiag(wi), sino(:), R, ig.mask);
		cost_lalm = pwls_cost(xlalm, A, Gdiag(wi), sino(:), R, ig.mask);
		cost_rlalm = pwls_cost(xrla, A, Gdiag(wi), sino(:), R, ig.mask);
	end
		fun = @(c) (c - cost_sqs(end)) / (cost_sqs(1) - cost_sqs(end));

	clf
	plot(	0:f.niter, fun(cost_sqs), '-s', ...
		0:f.niter, fun(cost_lalm), '-+', ...
		0:f.niter, fun(cost_rlalm), '-o')
	legend('SQS', 'LALM', 'Relaxed')
end
