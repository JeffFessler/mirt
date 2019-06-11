% fig_osem_mismatch.m
% examine effects of model mismatch on ML-OSEM
% Copyright 2002-6-10, Jeff Fessler, The University of Michigan

disp 'todo: NOT DONE!'

%
% generate data
%
if ~isvar('xfbp'), printm 'data'
	x = read_zubal_emis;
	[nx ny] = size(x);
	nb = 100;
	na = 10;

	arg = arg_pair('system', 12, 'nx', nx, 'ny', ny, ...
		'nb', nb, 'na', na, 'support', 'ellipse 0 0 60 60', ...
		'orbit', 180, 'orbit_start', 0, ...
		'pixel_size', 3, 'ray_spacing', 3, ...
		'obj2det_x', 2*nx, 'obj2det_y', 2*nx, ...
		'fwhm_detector', 0, ...
		'fwhm_collimator0', 0, ...
		'fwhm_slope', 0, ...
		'scale', 0);
prompt
end


%
% block object for block iterative
%
if ~isvar('Gb'), printm 'Gb'
	f.nblock = 5;
	Gb = Gblock(G, f.nblock);
prompt
end


%
% matlab iterations
%
if ~isvar('xmat'), printm 'matlab E-ML-OSEM'
	f.niter = 9;
	f.pixmax = 6;

	xinit = ones(size(xtrue));
	xmat = eml_osem(xinit(mask(:)), Gb, yi, ci, ri, f.niter, f.pixmax, ...
		[], 1.0);	% default fast precon matches aspire
	xmat = embed(xmat, mask);
	im clf, im(121, xmat, 'Matlab E-ML-OSEM iterations'), cbar horiz
prompt
end

if ~has_aspire, return, end

%
% aspire iterations
%
if ~isvar('xasp'), printm 'aspire E-ML-OSEM'

	f.init	= [f.dir 'init.fld'];
	f.out	= [f.dir 'out.fld'];
	fld_write(f.init, xinit, 'check', 0)
	if exist(f.out, 'file'), delete(f.out), end

	if ~isvar('f.sys_type')
		f.sys_type = ['2z@' f.wtr '@-'];
	end
	if ~isvar('n.n2')
		n.n2 = n.a;	% kludge
	end
	f.saver	= 'stack,1';
	f.alg	= sprintf('osemc,fast,%d,%d,1', f.nblock, n.n2);
	f.penal	= '-';
	f.method = sprintf('@%d@%s@%s', f.niter-1, f.alg, f.penal);
	f.com = sprintf(['i -chat 0 empl3 %s %s  %s %s 1 %s 1 %s -' ...
			' %s %s 0 1 %g 0 -'], ...
		f.out, f.init, f.yi, f.ci, f.ri, f.sys_type, ...
		f.method, f.saver, f.pixmax);
	os_run(f.com)

	xasp = double(fld_read(f.out));
	im(122, xasp, 'Aspire E-ML-OSEM iterations'), cbar horiz
end


if 1
	t = vcorrcoef(xasp, xmat);
	printf('corr. %g,%g', t, t-1)

	im clf, im(221, xmat, 'xhat matlab')
	im(222, xasp, 'xhat aspire')
	im(223, (xasp-xmat)/max(xmat(:)), 'aspire-matlab'), cbar

	t1 = eql_obj(xmat, G, yi(:), ci(:), ri(:), [], mask);
	t2 = eql_obj(xasp, G, yi(:), ci(:), ri(:), [], mask);

	if im
		subplot(224)
		plot(0:f.niter-1, t1-t1(1), '-o', 0:f.niter-1, t2-t1(1), '-x')
		xlabel iteration, ylabel '\Phi change', legend('mat', 'asp', 4)
		title(sprintf('E-ML-OSPS, Nsubset=%d', f.nblock))
	end
end
