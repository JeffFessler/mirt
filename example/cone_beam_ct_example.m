% cone_beam_ct_example.m
% Illustrate cone-beam X-ray CT image reconstruction via FDK and iterative.
% Uses a tiny system size because cone-beam projection can be slow.
% Copyright 2005-1-21, Jeff Fessler, University of Michigan

%% First run the FDK example.  It generates the true image xtrue
% and noiseless projection views "proj" and noisy data "yi"
% and generates (noisy) FDK recon "xfdk" for comparison / initialization.
if ~isvar('xfdk')
	bi = 1e6; % 1M photons / ray
	ri = 0; % no scatter etc. for now
	dfs = inf; % flat!
	feldkamp_example
prompt
end



%% generate system matrix needed for iterative reconstruction
if ~isvar('A'), printm 'A'
	f.nz_pad = 18; % add this many slices to each side.  todo: how many?
	ig_pad = ig.expand_nz(f.nz_pad);
	A = Gcone(cg, ig_pad, 'type', 'sf2', 'class', 'Fatrix');
%	f.sys_type = aspire_pair(cg, ig_pad, 'system', '3l'); % old way
%	A = Gtomo3(f.sys_type, ig_pad.mask, ig.nx, ig.ny, ig_pad.nz, ...
%		'chat', 0, 'permute213', true, 'checkmask', im&0);
end



% block object for ordered-subsets iterations
if ~isvar('Ab'), printm 'Ab'
	f.nblock = 8;
	Ab = Gblock(A, f.nblock);
end


if 1 % padding functions to help with long object problem
	f.rep = @(x) cat(3, repmat(x(:,:,1), [1 1 f.nz_pad]), x, ...
		repmat(x(:,:,end), [1 1 f.nz_pad]));
	f.pad = @(x) cat(3, zeros(ig.nx, ig.ny, f.nz_pad), x, ...
		zeros(ig.nx, ig.ny, f.nz_pad));
end


% check discrete vs analytical projection (note long object problem)
if 0, printm 'proj check'
	for ii=1:1
		cpu etic
		pp = Ab * f.rep(xtrue);
		cpu etoc 'proj time'
	end
	nrms(pp, proj)
	im clf, im_toggle(proj(:,:,1:12:end), pp(:,:,1:12:end), [0 4.4])
prompt
end


%% regularization object
if ~isvar('R'), printm 'regularizer'
	f.l2b = 2^4.5;
	f.delta = 100/1000;
	R = Reg1(ig_pad.mask, 'type_denom', 'matlab', ...
                'pot_arg', {'hyper3', f.delta}, 'beta', 2^f.l2b);
	W = Gdiag(yi);
	if 1 % check spatial resolution (away from edges)
		psf = qpwls_psf(A, R, 1, ig_pad.mask, W, 'fwhmtype', 'profile');
	prompt
	end
end


if ~isvar('si'), printm 'log sinogram'
	si = log(bi ./ yi); % log sinogram
	xinit = ig_pad.maskit(f.rep(xfdk)); % initalize iterations
	im(ig_pad.embed(xinit))
end


if 0 % todo: temporary lines for sam
	ytest = A * xtrue;
	xtest = A' * ytest
	wi = W.arg.diag % W = diag(wi)
	save sam1.mat xtest ytest xtrue cg ig si wi xinit
end


%% OS-SQS iterations for PWLS
if ~isvar('xos'), printm 'start os-sqs pwls iterations'
	f.niter_os = 10;
	xos = pwls_sqs_os(xinit, Ab, reshaper(si, '2d'), R, ...
			'wi', reshaper(W.arg.diag, '2d'), 'niter', f.niter_os);
	xos = ig_pad.embed(xos);
	im(xos), cbar
end


%% CG iterations
if ~isvar('xpcg1'), printm 'PWLS reconstruction with CG (no preconditioner)'
	f.niter1 = 20;
	[xpcg1 cost_pcg1] = pwls_pcg1(xinit, A, W, si(:), R, ...
		'userfun', @userfun_cost1, ...
		'niter', f.niter1, 'stop_threshold', 1e-3);
	xpcg1 = ig_pad.embed(xpcg1);
	im(xpcg1), cbar
prompt
end


% circulant preconditioner based on middle voxel and quadratic regularizer
if ~isvar('pre2'), printm 'circulant preconditioner'
	Rq = Reg1(ig_pad.mask, 'beta', 2^f.l2b); % quadratic regularizer
	pre2 = qpwls_precon('circ0', {A, W}, Rq.C, ig_pad.mask); % preconditioner
end


%% PCG
if 0 || ~isvar('xpcg2'), printm 'PWLS with PCG and circulant preconditioner'
	[xpcg2 cost_pcg2] = pwls_pcg1(xinit, A, W, si(:), R, ...
		'userfun', @userfun_cost1, ...
		'niter', f.niter1, 'stop_threshold', 1e-3, 'precon', pre2);
	xpcg2 = ig_pad.embed(xpcg2);
	im(xpcg2), cbar
prompt
end


if 0 % compare "convergence rate" (of cost) without and with preconditioning
	if im
		clf
		plot(1:f.niter1, cost_pcg1, '-o', 1:f.niter1, cost_pcg2, '-x')
		xlabel 'iteration', ylabel 'cost'
		legend('CG (no precon)', 'PCG2 (circulant precon)')
	end
return
end


% reshape data to be "2d arrays" for OS iterations (subset over last dim)
if ~isvar('os_data'), printm 'os_data'
	if isscalar(bi) && isscalar(ri)
		os_data = reshaper(yi, '2d');
		os_data = {os_data, ...
			bi * ones(size(os_data)), ri * ones(size(os_data))};
	else
		os_data = {reshaper(yi, '2d'), reshaper(bi, '2d'), ...
			reshaper(ri, '2d')}; % all data as 2d arrays
	end
end


%% OS-SPS iterations for transmission penalized likelihood
if ~isvar('xpl'), printm 'start tpl os-sqs iterations'
	f.niter = 20;
	xs = tpl_os_sps(xinit, Ab, os_data{:}, R, 1+f.niter);
	xs = ig_pad.embed(xs);
	xpl = xs(:,:,:,end);
	im(xpl)
end


%% compare FDK vs iterative results
if 1
	iz_good = f.nz_pad + [1:ig.nz]; % original slices
	xpcg1_good = xpcg1(:,:,iz_good);
	xpl_good = xpl(:,:,iz_good);
	xos_good = xos(:,:,iz_good);
	show = @(x) xlabelf('rms = %.1f HU', rms(col(xtrue-x))/0.01*1000);
	clim = [0 0.02];
	im plc 2 3
	im(1, xtrue, 'true', clim), cbar
	im(2, xfdk, 'FDK', clim), cbar
	show(xfdk)
	im(4, xpcg1_good, 'PWLS-PCG', clim), cbar
	show(xpcg1_good)
	im(5, xpl_good, 'PL-OS', clim), cbar
	show(xpl_good)
	im(6, xos_good, 'PWLS-OS', clim), cbar
	show(xos_good)
end
