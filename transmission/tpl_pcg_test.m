% tpl_pcg_test.m
% compare T-PL-PCG vs T-PL-OS-SPS
% Copyright 2004-2-1, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'), printm 'setup tpl_os_sps_test'
	tr_test_setup
prompt
end

if ~isvar('R'), printm 'R'
	f.l2b = 18;
	f.delta = 0.001;
%	f.pot = 'quad';
	f.pot = 'huber';
	R = Robject(ig.mask, 'type_denom', 'matlab', ...
		'potential', f.pot, 'beta', 2^f.l2b, 'delta', f.delta);
prompt
end


%
% pcg
%
if ~isvar('xpcg'), printm 'pcg'
	f.niter = 24+1;
	xinit = xfbp;

	M = 1; % no preconditioner
	if 0 % diagonal preconditioner
		ej = ig.unitv(floor(ig.nx/2)+1, floor(ig.ny/2)+1);
		gij_max = max(col(G*ej))
		curvi = trl_curvature(yi, bi, ri, 0, 'pc');
		M = 1 ./ (gij_max * (G' * curvi(:)) + R.diag(R));
		M = spdiag(M);
	end
	data = {yi(:), bi(:), ri(:)};
	xpcg = pl_pcg_qs_ls(xinit(ig.mask), G, data, ...
		@trl_dercurv, R, 'precon', M, 'niter', f.niter, ...
		'curvtype', 'oc', 'isave', 'all');
	xpcg = ig.embed(xpcg);
	im clf, im(xpcg, 'pcg iterations'), cbar
prompt
end


%
% os-sps
%
if ~isvar('Gb'), printm 'Gb'
	f.nblock = 5;
	Gb = Gblock(G, f.nblock);
prompt
end


if ~isvar('xos'), printm 'os'
	xos = tpl_os_sps(xinit(ig.mask), Gb, yi, bi, ri, R, ...
			f.niter, [-inf inf], 'pc');
	xos = ig.embed(xos);
	im clf, im(xos, 'OS'), cbar
prompt
end

if 1
	ipcg = xpcg(:,:,end);
	ios = xos(:,:,end);
	t = vcorrcoef(ipcg, ios);
	printf('PCG vs OS corr. %g,%g', t, t-1)
	printf('PCG vs OS max %% diff = %g', max_percent_diff(ipcg, ios))

	im clf, im(221, ios, 'xos'), cbar
	im(222, ipcg, 'xpcg'), cbar
	im(223, 100*(ipcg-ios)/max(col(ios)), '100(pcg-os)/max(os)'), cbar

	t1 = tpl_obj(xpcg, G, yi(:), bi(:), ri(:), R, ig.mask);
	t2 = tpl_obj(xos, G, yi(:), bi(:), ri(:), R, ig.mask);

	if im
		subplot(224)
		i0 = 80:f.niter-1;
		i0 = 0:f.niter-1;
		plot(i0, t1(i0+1)-t1(1), '-o', i0, t2(i0+1)-t1(1), '-x')
		xlabel iteration, ylabel '\Psi change', legend('PCG', 'OS', 4)
%		title(sprintf('T-PL-OSPS, Nsubset=%d', f.nblock))
	end
end
