% eml_em_test.m
% compare aspire and matlab E-ML-EM
% Copyright Jan 1998	Jeff Fessler, The University of Michigan


%
% generate data
%
if ~isvar('yi'), printm 'setup for eml_em_test'
	if has_aspire
		f.dir	= test_dir;
		f.wtf	= [f.dir 't,g.wtf'];
		f.yi	= [f.dir 'yi.fld'];
		f.ci	= [f.dir 'ci.fld'];
		f.ri	= [f.dir 'ri.fld'];
		f.mask	= [f.dir 'mask.fld'];
	end
	em_test_setup
prompt
end

%
% matlab iterations
%
if ~isvar('xmat'), printm 'matlab E-ML-EM'
	f.niter = 9;
	xinit = ig.ones; % uniform
	xmat = eml_em(xinit(ig.mask), G, yi(:), ci(:), ri(:), [], f.niter);
	xmat = ig.embed(xmat);

	im clf, im(xmat, 'matlab E-ML-EM iterates')
prompt
end

if ~has_aspire, return, end

%
% aspire iterations - this is for Fessler's testing only!
%
if ~isvar('xasp'), printm 'aspire E-ML-EM'

	f.init	= [f.dir 't,init.fld'];
	f.out = [f.dir 't,out.fld'];
	fld_write(f.init, xinit, 'check', 0)
	if (exist(f.out) == 2), delete(f.out), end

	f.alg = 'em,1';
	f.saver = '-';
	f.saver = 'stack,1';
	f.penal = '-';
	f.method = sprintf('@%d@%s@%s', f.niter-1, f.alg, f.penal);
	f.scaleinit = 0;

	if 1
		f.com = sprintf('i -chat 0 empl2 %s %s  %s %s %s 1 %s %s  %s %s 0 1e9 %d -', ...
			f.out, f.init, f.yi, f.ci, f.ri, f.wtf, ...
			f.mask, f.method, f.saver, f.scaleinit);
		os_run(f.com)

		xasp = double(fld_read(f.out));
	end
end

if 1
	im(221, xmat, 'xhat matlab')
	im(222, xasp, 'xhat aspire')
	im(223, (xasp-xmat)/max(xmat(:)), 'aspire-matlab')

	t = vcorrcoef(xasp, xmat);
	printf('corr. %g,%g', t, t-1)

	t1 = eql_obj(xmat, G, yi(:), ci(:), ri(:), [], ig.mask);
	t2 = eql_obj(xasp, G, yi(:), ci(:), ri(:), [], ig.mask);

	if im
		subplot(224)
		plot(0:f.niter-1, t1-t1(1), '-o', 0:f.niter-1, t2-t1(1), '-x')
		legend('mat', 'asp', 4), xlabel iteration, ylabel objective
	end
end
