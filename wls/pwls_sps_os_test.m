% pwls_sps_os_test.m
% Test 2d PWLS-SPS-OS with nonquadratic penalty
%
% Copyright 2002-2-13, Jeff Fessler, University of Michigan

% generate data
if ~isvar('yi'), printm 'generating data for pwls_sps_os_test'
	f.dir	= test_dir;
	f.wtf	= [f.dir 't.wtf'];
	f.yi	= [f.dir 'yi.fld'];
	f.wi	= [f.dir 'wi.fld'];
	f.wtr = strrep(f.wtf, '.wtf', ',row.wtf');
	ig = image_geom('nx', 32, 'ny', 30, 'fov', 500);
	tmp = ig.circ > 0;
	tmp(:,[2 end-1]) = 0; % double border of zeros
	ig.mask = tmp;
%	clf, im(ig.mask), return

	if 0 % tiny case so exact solution for quadratic
		ig = image_geom('nx', 12, 'ny', 14, 'fov', 500);
%		n.x=12; n.y=14; n.a=30; n.b=24;
		f.ray_spacing = 0.5;
%		f.ray_pix = 0.5;
	end

	em_wls_test_setup

	W = diag_sp(wi(:));
	if has_aspire
		Awtr = Gtomo2_wtmex(f.wtr); % need row grouped!
	else
		Awtr = G; % Gtomo2_strip
	end
	clear G
prompt
end


% block system object
if ~isvar('Ab'), printm 'make block system'
	f.nblock = 8;
%	f.nblock = 1; % for testing vs true xhat
	Ab = Gblock(Awtr, f.nblock);
prompt
end


% regularizer
if ~isvar('R'), printm 'make R'
%	kappa = ig.embed(sqrt((A' * wi(:)) ./ sum(A)'));
	kappa = ig.mask; % to match plain aspire penalty
	im(8, kappa, 'kappa'), cbar

	f.l2b = 0;
	f.nbrs = 8;
	f.delta = 0.5;
	f.pot = 'huber';
	if 0, printm 'testing quadratic case'
		f.pot = 'quad';
		f.delta = [];
	end

	if 0
		R = Robject(ig.mask, 'edge_type', 'leak', 'beta', 2^f.l2b, ...
			'potential', f.pot, 'delta', f.delta, 'type_denom', 'matlab');
		im(9, ig.shape(R.wt), 'R.wt'), cbar
	else
		R = Reg1(ig.mask, 'edge_type', 'leak', 'beta', 2^f.l2b, ...
			'pot_arg', {f.pot, f.delta}, 'type_denom', 'matlab', ...
			'type_penal', 'mat'); % to match aspire
		im subplot 9
		R.wt.show; cbar
	end

prompt
	clear xhat xmat xasp
end


% exact solution
if ~isvar('xhat') && streq(f.pot, 'quad'), printm 'exact xhat'
	Atmp = Awtr(:,:);
	F = full(Atmp' * W * Atmp);
	C1 = R.C1(:,:);
	Rh = full(C1' * spdiag(R.wt) * C1);	% penalty Hessian
	if 0 % verify Hessian matches quadratic penalty
		x = xtrue(ig.mask);
		[R.penal(R, x) 0.5 * x' * Rh * x]
	return
	end
	printm('cond F = %g', cond(F))
	printm('cond(F+R) = %g', cond(F+Rh))
	xhat = (F + Rh) \ (Atmp' * W * yi(:));
	xhat = ig.embed(xhat);
prompt
	clear F Rh Atmp C1
end


% matlab iterations
if ~isvar('xmat'), printm 'matlab PWLS-SPS-OS'

	xinit = max(xfbp,0); % FBP initial image

	if streq(f.pot, 'quad')
		f.pixmin = -inf;
		f.niter = 25+1;
	else
		f.pixmin = 0;
		f.niter = 8+1;
	end
	f.pixmax = 1e9;

	denom = []; gi = [];
	if 0
		xmat = pwls_sps_os(xinit(ig.mask), yi, wi, ...
			Ab, R, f.niter, [f.pixmin f.pixmax], denom, gi, 1, 0);
	else
		xmat = pwls_sqs_os(xinit(ig.mask), Ab, yi, R, 'wi', wi, ...
			'niter', f.niter-1, 'pixmax', [f.pixmin f.pixmax], ...
			'denom', denom, 'aai', gi, 'isave', 'all');
	end
	xmat = ig.embed(xmat);

	im clf, im(xmat, 'matlab PWLS-SPS-OS')

	% examine matlab convergence to xhat
	if isvar('xhat'), printm 'matlab vs xhat'
		prompt
		err = reshape(xmat, ig.nx*ig.ny, f.niter) - repmat(xhat(:),1,f.niter);
		err = sqrt(sum(err.^2)) / norm(xhat(:));
		semilogy(err, '.-'), axis tight, title 'xmat-xhat'
	end

prompt
end

if ~has_aspire, return, end

% aspire iterations
if ~isvar('xasp'), printm 'aspire PWLS-SPS-OS'

	f.init	= [f.dir 'init.fld'];
	fld_write(f.init, xinit, 'check', 0)

	f.out	= [f.dir 'out.fld'];
	if exist(f.out, 'file'), delete(f.out), end

	if isinf(f.pixmin)
		f.nonneg = false;
	else
		f.nonneg = true;
	end
	if streq(f.pot, 'quad')
		f.penal	= sprintf('%g,quad,%d,-', f.l2b, f.nbrs/4);
	elseif streq(f.pot, 'huber')
		f.penal	= sprintf('%g,huber,%d,-,%g,ih,1', ...
			f.l2b, f.nbrs/4, f.delta);
	else
		error 'bad potential'
	end
	f.saver	= 'stack,1';
	f.alg	= sprintf('ospsc,%d,%d,1,0', f.nblock, sg. na);
	f.method = sprintf('@%d@%s@%s', f.niter-1, f.alg, f.penal);
	f.fitype = ['2z@' f.wtr '@-'];
	f.com = sprintf(['i -chat 0 pwls3 %s %s %s %s %s' ...
			' - %s %s 1 %d %g 0 -'], ...
		f.out, f.init, f.yi, f.fitype, f.wi, ...
		f.method, f.saver, f.nonneg, f.pixmax);

	os_run(f.com)

	xasp = fld_read(f.out);
	im clf, im(xasp, 'xasp')

	% examine aspire convergence to xhat
	if isvar('xhat'), printm 'aspire vs xhat'
		prompt
		err = reshape(xasp, ig.nx*ig.ny, f.niter) - repmat(xhat(:),1,f.niter);
		err = sqrt(sum(err.^2)) / norm(xhat(:));
		semilogy(err, '.-'), axis tight, title 'xasp-xhat'
	end
prompt
end

if 1
	equivs(xasp, xmat)

	im plc 2 2
	im(1, xmat, 'xhat matlab'), cbar
	im(2, xasp, 'xhat aspire'), cbar
	im(3, xasp-xmat, 'aspire-matlab'), cbar

	t0 = pwls_cost(xinit, Ab, W, yi(:), R, ig.mask);
	t1 = pwls_cost(xmat, Ab, W, yi(:), R, ig.mask);
	t2 = pwls_cost(xasp, Ab, W, yi(:), R, ig.mask);

	if im
		im subplot 4
		plot(0:f.niter-1, t1-t0, '-o', 0:f.niter-1, t2-t0, '-x')
		xlabel iteration, ylabel '\Phi change'
		ir_legend({'mat', 'asp'})
		titlef('PWLS-SPS-OS, Nsubset=%d', f.nblock)
	end
end
