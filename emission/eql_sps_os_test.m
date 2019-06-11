% eql_sps_os_test.m
% compare aspire and matlab versions of E-QPL-SPS-OS
% Copyright Apr 2000, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'), printm 'data'
	if has_aspire
		f.dir	= test_dir;
		f.wtf	= [f.dir 't,g.wtf'];
		f.wtr	= strrep(f.wtf, 'wtf', 'wtr');
		f.yi	= [f.dir 'yi.fld'];
		f.ci	= [f.dir 'ci.fld'];
		f.ri	= [f.dir 'ri.fld'];
		f.mask	= [f.dir 'mask.fld'];
	end
	em_test_setup
prompt
end

if ~isvar('Gb'), printm 'Gb'
	f.nblock = 5;
	Gb = Gblock(G, f.nblock);
prompt
end

if ~isvar('R'), printm 'R'
	f.nbrs = 8;
	f.l2b = 2;
	R = Robject(ig.mask, 'edge_type', 'leak', 'type_denom', 'aspire', ...
		'beta', 2^f.l2b);
prompt
end

%
% matlab iterations
%
if ~isvar('xmat'), printm 'matlab E-QPL-SPS-OS'
	f.niter = 9;
	f.pixmax = 6;
	curvs = {'pc', 'oc'};
	xinit = max(xfbp,0);
	for ic = 1:length(curvs)
		tmp = eql_sps_os(xinit(ig.mask), Gb, yi, ci, ri, R, ...
			f.niter, f.pixmax, curvs{ic});
		xmat{ic} = ig.embed(tmp);
	end
	im clf, im(xmat{1}, 'Matlab E-QPL-SPS-OS iterations')
prompt
end


if ~has_aspire, return, end

%
% aspire iterations
%
if ~isvar('xasp'), printm 'aspire E-QPL-SPS-OS'

	f.init	= [f.dir 'init.fld'];
	f.out	= [f.dir 'out.fld'];
	fld_write(f.init, xinit, 'check', 0)

	f.saver	= 'stack,1';
	f.fitype = ['2z@' f.wtr '@-'];

	for ic = 1:length(curvs)
		if exist(f.out, 'file'), delete(f.out), end
		f.alg = sprintf('ospsc,%s,%d,%d,1,0', ...
			curvs{ic}, f.nblock, sg.na);
		f.penal	= sprintf('%g,quad,%d,-', f.l2b, f.nbrs/4);
		f.method = sprintf('@%d@%s@%s', f.niter-1, f.alg, f.penal);
		f.com = sprintf(['i -chat 0 empl3 %s %s  %s %s 1 %s 1 %s -' ...
				' %s %s 0 1 %g 0 -'], ...
			f.out, f.init, f.yi, f.ci, f.ri, f.fitype, ...
			f.method, f.saver, f.pixmax);

		os_run(f.com)
		xasp{ic} = double(fld_read(f.out));
	end
end

if 1
	for ic = 1:length(xasp)
		t = vcorrcoef(xasp{ic}, xmat{ic});
		printf('corr. %g,%g', t, t-1)
		err = (xasp{ic} - xmat{ic}) / max(col(xmat{ic}));
		printf('Normalized error range %g %g', min(err(:)), max(err(:))),
	end

	im clf, im(221, xmat{1}, 'xhat matlab'), cbar
	im(222, xasp{1}, 'xhat aspire'), cbar
	im(223, (xasp{1}-xmat{1})/max(col(xmat{1})), 'aspire-matlab'), cbar

	t1 = eql_obj(xmat{1}, G, yi(:), ci(:), ri(:), R, ig.mask);
	t2 = eql_obj(xasp{1}, G, yi(:), ci(:), ri(:), R, ig.mask);

	if im
		subplot(224)
		plot(0:f.niter-1, t1-t1(1), '-o', 0:f.niter-1, t2-t1(1), '-x')
		xlabel iteration, ylabel '\Phi change', legend('mat', 'asp', 4)
		title(sprintf('E-QPL-SPS-OS, Nsubset=%d', f.nblock))
	end
end
