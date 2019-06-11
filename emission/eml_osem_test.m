% eml_osem_test.m
% compare aspire and matlab versions of E-ML-OSEM
% Copyright 2001-07-23, Jeff Fessler, The University of Michigan

% generate data
if ~isvar('xfbp'), printm 'data'
	if has_aspire
		f.dir	= test_dir;
		f.wtf	= [f.dir 't,g.wtf'];
		f.wtr	= strrep(f.wtf, 'wtf', 'wtr');
		f.yi	= [f.dir 'yi.fld'];
		f.ci	= [f.dir 'ci.fld'];
		f.ri	= [f.dir 'ri.fld'];
		f.mask	= [f.dir 'mask.fld'];
	end
	em_test_setup; f.is3b = 0;
%	em3_test_setup; f.is3b = 1;
prompt
end


% block object for block iterative
if ~isvar('Gb'), printm 'Gb'
	f.nblock = 5;
	Gb = Gblock(G, f.nblock);
end


% matlab iterations
if ~isvar('xmat'), printm 'matlab E-ML-OSEM'
	f.niter = 8;
	f.pixmax = 6;
	f.precon = 'fast';
	f.precon = 'classic';

	if f.is3b % fix: needs work!
		shaper = @(proj) reshape(proj, n.n1*n.n2, n.a);
	else
		shaper = @(x) x;
	end
	xinit = ig.ones;
	[xmat precon] = eml_osem(xinit(ig.mask), Gb, ...
			shaper(yi), shaper(ci), shaper(ri), ...
			'niter', f.niter, 'pixmax', f.pixmax, ...
			'precon', f.precon);
	precon = ig.embed(precon);
	xmat = ig.embed(xmat);
	im plc 1 2
	im(1, xmat, 'Matlab E-ML-OSEM iterations'), cbar horiz
prompt
end

if ~has_aspire, return, end

% aspire iterations
if ~isvar('xasp'), printm 'aspire E-ML-OSEM'

	f.init	= [f.dir 'init.fld'];
	f.out	= [f.dir 'out.fld'];
	fld_write(f.init, xinit, 'check', 0)
	if exist(f.out, 'file'), delete(f.out), end

	if ~isvar('f.sys_type')
		f.sys_type = ['2z@' f.wtr '@-'];
	end
	f.saver	= 'stack,1';
	f.alg	= sprintf('osemc,%s,%d,%d,1', f.precon, f.nblock, sg.na);
	f.penal	= '-';
	f.method = sprintf('@%d@%s@%s', f.niter, f.alg, f.penal);
	f.com = sprintf(['i -chat 5 empl3 %s %s  %s %s 1 %s 1 %s %s' ...
			' %s %s 0 1 %g 0 -'], ...
		f.out, f.init, f.yi, f.ci, f.ri, f.sys_type, f.mask, ...
		f.method, f.saver, f.pixmax);
	os_run(f.com)

	xasp = double(fld_read(f.out));
	im(2, xasp, 'Aspire E-ML-OSEM iterations'), cbar horiz
end

if 0 % test aspire preconditioner
	pasp = squeeze(fld_read('p.fld'));
	im plc 1 2
	im(1, precon, 'precon'), cbar
	im(2, pasp, 'aspire'), cbar
	max_percent_diff(precon, pasp)
	minmax(precon-pasp)
return
end


if 1
	t = vcorrcoef(xasp, xmat);
	printf('corr. %g,%g', t, t-1)

	im plc 2 2
	im(1, xmat, 'xhat matlab'), cbar
	im(2, xasp, 'xhat aspire'), cbar
	im(3, (xasp-xmat)/max(xmat(:)), 'aspire-matlab'), cbar

	t1 = eql_obj(xmat, G, yi(:), ci(:), ri(:), [], ig.mask);
	t2 = eql_obj(xasp, G, yi(:), ci(:), ri(:), [], ig.mask);

	if im
		subplot(224)
		plot(0:f.niter, t1-t1(1), '-o', 0:f.niter, t2-t1(1), '-x')
		xlabel iteration, ylabel '\Phi change', legend('mat', 'asp', 4)
		title(sprintf('E-ML-OSEM, Nsubset=%d', f.nblock))
	end
end
