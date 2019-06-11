% pwls_gca_test.m
% Test 2d penalized weighted least-squares with nonquadratic penalty
%
% Copyright Mar 1999, Jeff Fessler, University of Michigan

% generate data
if ~isvar('xtrue'), printm 'xtrue'
	ig = image_geom('nx', 12, 'ny', 10, 'dx', 1);
	ig.mask = ig.circ > 0;
	xtrue = ig.circ(2);
	im plc 3 2
	im(1, xtrue, 'xtrue')
	im(2, ig.mask, 'mask')
prompt
end


% form A matrix
if ~isvar('A'), printm 'A'
	sg = sino_geom('par', 'nb', 14, 'na', 8, 'dr', 0.5);
	rng(0)
	A = Gtomo2_strip(sg, ig);
	A = A(:,:); % need matrix for kappa later
%	A = A(:,ig.mask(:));
%	A = full(A); A = sparse(A(:,mask));

	if has_mex_jf
		dir = test_dir;
		f.wtf = [dir 't2,g.wtf'];
		if exist(f.wtf), delete(f.wtf), end
		Atmp = zeros(sg.nb * sg.na, ig.nx * ig.ny);
		Atmp(:,ig.mask) = A;
		wtf_write(f.wtf, Atmp, ig.nx, ig.ny, sg.nb, sg.na);
		clear Atmp
	end

	if im, clf, spy(A(:,:)), title('A'), end
%	printm('condest(At*A) = %g', condest(A'*A))
prompt
end


if ~isvar('yi'), printm 'yi'
	proj = sg.shape(A * xtrue(ig.mask));
	count = 1e4;
	ci = count / sum(proj(:)) * sg.ones;
	ytrue = ci .* proj;
	randpercent = 10;
	ri = randpercent / 100 * mean(ytrue(:));
	rng(0)
	yi = poisson(ytrue + ri);

	im(3, ytrue, 'ytrue')
	im(4, yi, 'yi')
prompt
end


% "exact" unconstrained *quadratic penalty* solution
if ~isvar('xhat'), printm 'xhat'
	pivot = max((yi - ri) ./ ci, 0);
	pivot = dsingle(pivot);
	im(1, proj, 'proj')
	im(3, pivot, 'pivot')
	im(5, pivot-proj, 'pivot-proj')

	ybarhat = ci .* pivot + ri; % positive because ri > 0
	nder1 = ci .* (1 - yi ./ ybarhat);
%	nder1 = zeros(n.dy);
	nder1 = dsingle(nder1);
	minmax(nder1)
	nder2 = yi .* (ci ./ ybarhat).^2;
%	nder2 = ones(n.dy);
	nder2 = dsingle(nder2);
	W = diag_sp(nder2(:));

	kappa = ig.embed(sqrt(((A.*A)' * nder2(:)) ./ sum(A.*A)'));
	im(2, kappa, 'kappa')

	% b2info
	f.l2b = -1; C = sqrt(2^f.l2b) * Csparse('b2info', kappa, ig.mask, 0);
	f.penal = sprintf('%g,quad,1,b2info', f.l2b);

	% standard
	f.l2b = 12; C = sqrt(2^f.l2b) * Csparse('maskleak', ig.mask, 0);
	f.penal = sprintf('%g,quad,1,-', f.l2b);

	% nonquadratic case
	if 1
		f.l2b = 6;
		C = Csparse('maskleak', ig.mask, 0);
		f.delta = 1.5;
		f.ptype = 'lange3';
		f.nsub = 3;
		[Rarg.pstring, Rarg.wstring, Rarg.cstring] = ...
			rp_string(f.ptype, 2^f.l2b, f.delta);
		f.penal = sprintf('%g,%s,1,-,%g,ih,%d', ...
			f.l2b, f.ptype, f.delta, f.nsub)
	end

	if 0
		t = reshape(C * xtrue(ig.mask), [n.dx 2]);
		im(6, t)
	return
	end

	F = double(full(A' * W * A));
	sprintf('cond = %g', cond(full(F + C'*C)))
	backs = A' * (W * pivot(:) - nder1(:));
	xhat = (F + C'*C) \ backs;
	xhat = ig.embed(xhat);
	backs = ig.embed(backs);
	im(4, backs, 'backs')
	im(6, xhat, 'xhat')
prompt
end

% do pcg to compare to xhat

if 1
%	f.backs	= [dir 't2,backs.fld';
	f.nder1	= [dir 't2,nder1.fld'];
	f.nder2	= [dir 't2,nder2.fld'];
	f.pivot	= [dir 't2,pivot.fld'];
	f.init	= [dir 't2,init.fld'];
	f.mask	= [dir 't2,mask.fld'];
end

if ~exist(f.mask, 'file') || 1
%	xinit = xtrue > 0;	%	f.init = '-'
%	xinit = xhat;
	xinit = xtrue;
	xinit = dsingle(xinit);
	delete([dir 't2,*.fld'])
	fld_write(f.init, xinit)
%	fld_write(f.backs, backs)
	fld_write(f.pivot, pivot)
	fld_write(f.nder1, nder1)
	fld_write(f.nder2, nder2)
	fld_write(f.mask, ig.mask)
end

if 0
	t = reshape(C * xinit(mask), [n.dx 2]);
	im(1, t), sum(t(:).^2/2)
prompt
end

	f.niter = 20;
	n.gx = 3;
	n.gy = 3;

% Run Matlab
if ~isvar('groups')
	groups = group2d(ig.nx, ig.ny, [n.gx n.gy], ig.mask, 1);
	groups = logical(groups);
end
if ~isvar('xmat')
	xmat = pwls_gca(A, W, double(pivot(:)), double(nder1(:)), ...
		double(xinit(ig.mask)), ...
		C, f.niter+1, f.nsub, groups, Rarg.wstring, ig.mask, 1);
	xmat = ig.embed(xmat);
	im clf, im(1, xmat)
prompt
end


% Run ASPIRE
if ~has_aspire, return, end

if 1
	f.out = [dir 't2,out.fld'];
	f.alg = 'ca,1.0,raster1';
	f.alg = 'cg,none';
	f.alg = sprintf('ga,1.0,%d,%d', n.gx, n.gy);
	f.saver = '-';
	f.saver = 'stack,1';
	f.method = sprintf('@%d@%s@%s', f.niter, f.alg, f.penal);
	f.scaleinit = 0;

	if 1
		delete(f.out)
		f.com = sprintf('i pwls2 %s %s  %s %s %s  %s %s  %s %s 1 0 1e9 %d -', ...
			f.out, f.init, f.pivot, f.nder1, f.nder2, f.wtf, ...
			f.mask, f.method, f.saver, f.scaleinit);
		os_run(f.com)

		xasp = fld_read(f.out);

		im(221, xmat, 'xhat matlab')
		im(222, xasp, 'xhat aspire')
		im(212, xasp(:,:,end)-xmat(:,:,end), 'aspire-matlab')
%	prompt
	end
end

if exist('pwls_obj.m', 'file') % jf only
	%cost = @(x) pwls_cost(x, A, W, yi, R, mask)
	cost = @(x) pwls_obj(double(x(ig.mask)), ...
		A, W, pivot(:), nder1(:), C, Rarg.pstring, 0, 1);

	pr 'cost(xinit)'
	pr 'cost(xmat(:,:,end))'
	pr 'cost(xasp(:,:,end))'
end

printf('norm. error: %g%%', norm(xasp(:)-xmat(:)) / norm(xmat(:)) * 100)
t = xasp(:)'*xmat(:) / norm(xasp(:)) / norm(xmat(:));
printm('corr. %g, %g', t, t-1)
