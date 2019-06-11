% eml_sps_os_test.m
% compare aspire and matlab versions of E-ML-SPS-OS
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
%	em3_test_setup, n.a = n.n2; % kludge
prompt
end

if ~isvar('Gb'), printm 'make Gb'
	f.nblock = 5;
	Gb{1} = Gblock(G, 1);
	Gb{2} = Gblock(G, f.nblock);
prompt
end

%
% matlab iterations
%
if ~isvar('xmat'), printm 'matlab E-ML-SPS-OS'
	f.niter = 9;
	f.pixmax = 6;

	f.curvs = {'oc', 'pc'};

	xinit = ones(size(xtrue));
	for ic = 1:length(f.curvs)
		tmp = eql_sps_os(xinit(ig.mask), Gb{ic}, yi, ci, ri, [], ...
			f.niter, f.pixmax, f.curvs{ic});
		xmat{ic} = ig.embed(tmp);
	end
	im clf, im(xmat{1}, 'Matlab E-ML-SPS-OS iterations')
prompt
end

if ~has_aspire, return, end

%
% aspire iterations
%
if ~isvar('xasp'), printm 'aspire E-ML-SPS-OS'

	f.init	= [f.dir 'init.fld'];
	f.out	= [f.dir 'out.fld'];
	fld_write(f.init, xinit, 'check', 0)

	f.saver	= 'stack,1';
	f.fitype = ['2z@' f.wtr '@-'];

	for ic = 1:length(f.curvs)
		if exist(f.out, 'file'), delete(f.out), end
		f.alg = sprintf('ospsc,%s,%d,%d,1,0', f.curvs{ic}, Gb{ic}.nblock, sg.na);
		f.penal	= sprintf('%g,quad,0,-', -100);
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
	for ic = 1:length(f.curvs)
		t = vcorrcoef(xasp{ic}, xmat{ic});
		printf('corr. [curv=%s block=%d] %g,%g', ...
			f.curvs{ic}, Gb{ic}.nblock, t, t-1)
	end
	ic = 2;

	im clf, im(221, xmat{ic}, 'xhat matlab')
	im(222, xasp{ic}, 'xhat aspire')
	im(223, (xasp{ic}-xmat{ic})/max(col(xmat{ic})), 'aspire-matlab'), cbar

	t1 = eql_obj(xmat{ic}, G, yi(:), ci(:), ri(:), [], ig.mask);
	t2 = eql_obj(xasp{ic}, G, yi(:), ci(:), ri(:), [], ig.mask);

	if im
		subplot(224)
		plot(0:f.niter-1, t1-t1(1), '-o', 0:f.niter-1, t2-t1(1), '-x')
		xlabel iteration, ylabel '\Phi change', legend('mat', 'asp', 4)
		title(sprintf('E-ML-SPS-OS, Nsubset=%d', Gb{ic}.nblock))
	end
end
