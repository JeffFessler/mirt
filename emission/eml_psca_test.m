% eml_psca_test.m
% compare aspire and matlab versions of E-ML-PSCA
% Copyright Apr 2000, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'),	disp 'setup eml_psca_test'
	if has_aspire
		f.dir	= test_dir;
		f.wtf	= [f.dir 't,g.wtf'];
		f.yi	= [f.dir 'yi.fld'];
		f.ci	= [f.dir 'ci.fld'];
		f.ri	= [f.dir 'ri.fld'];
		f.mask	= [f.dir 'mask.fld'];
	end
	ig = image_geom('nx', 16, 'ny', 14, 'fov', 500);
	em_test_setup
prompt
end


%
% matlab iterations
%
if ~isvar('xmat'),	disp 'matlab E-ML-PSCA'
	f.niter = 9;

	xinit = ig.ones;
	x = xinit(ig.mask);

	xmat = zeros(sum(ig.mask(:)), f.niter);
	xmat(:,1) = x;
	for ii=2:f.niter
		x = eml_psca(x, G, yi, ci, ri);
		xmat(:,ii) = x;
	end
	xmat = ig.embed(xmat);

	im clf, im(xmat)
prompt
end

if ~has_aspire, return, end

%
% aspire iterations
%
if ~isvar('xasp'),	disp 'aspire E-ML-PSCA'

	f.init	= [f.dir 't,init.fld'];
	fld_write(f.init, xinit, 'check', 0)
	f.out = [f.dir 't,out.fld'];
	if (exist(f.out) == 2), delete(f.out), end

	f.scaleinit = 0;

	f.saver = 'stack,1';
	f.alg = 'psca,1,raster1';
	f.penal = sprintf('%g,quad,0,-', -100);
	f.method = sprintf('@%d@%s@%s', f.niter-1, f.alg, f.penal);
%	f.fitype = ['2z@' f.wtr '@-'];
%	f.com = sprintf(['i -chat 0 empl3 %s %s  %s %s 1 %s 1 %s -' ...
%			' %s %s 1 1 1e30 0 -'], ...
%		f.out, f.init, f.yi, f.ci, f.ri, f.fitype, ...
%		f.method, f.saver);

	f.com = sprintf(['i -chat 0 empl2 %s %s  %s %s %s 1 %s -' ...
			' %s %s 0 1e30 0 -'], ...
		f.out, f.init, f.yi, f.ci, f.ri, f.wtf, ...
		f.method, f.saver);

		os_run(f.com)
		xasp = double(fld_read(f.out));
end

	t = vcorrcoef(xasp, xmat);
	printf('corr. %g,%g', t, t-1)

	t1 = eql_obj(xmat, G, yi(:), ci(:), ri(:), [], ig.mask);
	t2 = eql_obj(xasp, G, yi(:), ci(:), ri(:), [], ig.mask);

if 1 && im
	im clf, im(221, xmat, 'xhat matlab'), cbar
	im(222, xasp, 'xhat aspire'), cbar
	im(223, (xasp-xmat)/max(xmat(:)), 'aspire-matlab'), cbar

	subplot(224)
	plot(0:f.niter-1, t1-t1(1), '-o', 0:f.niter-1, t2-t1(1), '-x')
	legend('mat', 'asp', 4)
prompt
end
