% eql_os_emdp_test.m
% compare aspire and matlab E-QPL-OS-EMDP (De Pierro)
% Copyright Nov 1998, Jeff Fessler, University of Michigan

%
% gen data
%
if ~isvar('yi'), printm 'data'
	if has_aspire
		f.dir	= test_dir;
		f.wtf   = [f.dir 't,g.wtf'];
		f.wtr	= strrep(f.wtf, 'wtf', 'wtr');
		f.yi	= [f.dir 'yi.fld'];
		f.ci	= [f.dir 'ci.fld'];
		f.ri	= [f.dir 'ri.fld'];
		f.mask	= [f.dir 'mask.fld'];
	end
	em_test_setup
end

if ~isvar('Ab'), printm 'Ab'
	f.nblock = 4;
	Ab = Gblock(G, f.nblock);
prompt
end

if ~isvar('R'), printm 'R'
	f.nbrs = 4; f.offsets = [1 ig.nx];
%	f.nbrs = 8; f.offsets = [1 ig.nx ig.nx-1 ig.nx+1]; % not done in aspire?
	f.l2b = 0;
	R = Robject(ig.mask, 'edge_type', 'leak', ... % to match aspire
		'offsets', f.offsets, ...
		'beta', 2^f.l2b, 'potential', 'quad', ...
		'type_denom', 'aspire');
prompt
end

%
% matlab iterations
%
if ~isvar('xmat'), printm 'matlab E-PL-EMDP'
	f.niter = 9;

	xinit = max(xfbp, 0.001*ig.mask); % requires x>0

	xmat = eql_os_emdp(xinit(ig.mask), Ab, yi, ci, ri, R, ...
		'niter', f.niter, 'hds', 1, 'chat', 0);
	xmat = ig.embed(xmat);

	im clf, im(xmat, 'Matlab EQL-OS-EMDP')
prompt
end

if ~has_aspire, return, end

%
% aspire iterations
%
if ~isvar('xasp'), printm 'aspire E-PL-EMDP'
	f.init	= [f.dir 't,init.fld'];
	f.out = [f.dir 't,out.fld'];
	fld_write(f.init, xinit, 'check', 0)
	if (exist(f.out) == 2), delete(f.out), end

	f.scaleinit = 0;

	if 1
		f.saver = 'stack,1';
		f.alg = sprintf('osdpc,%d', f.nblock);
		f.penal = sprintf('%g,quad,%d,-', f.l2b, f.nbrs/4);
		f.method = sprintf('@%d@%s@%s', f.niter, f.alg, f.penal);
		f.com = sprintf(['i -chat 0 empl2 %s %s  %s %s %s 1 %s -' ...
			' %s %s 0 1e30 0 -'], ...
			f.out, f.init, f.yi, f.ci, f.ri, f.wtr, ...
			f.method, f.saver);
	else
		f.out = 't,out,%02d.fld';
		f.alg = 'map-adp';
		f.saver = '-1';
		f.penal = '-';
		f.com = sprintf(['i -chat 0 em-2d %s %s  %s %s %s 1 %s ' ...
			' %s 0 %d %s %s  1 1 1 -'], ...
			f.out, f.init, f.yi, f.ci, f.ri, f.wtf, ...
			f.alg, f.niter, f.saver, f.penal);
	end

	os_run(f.com)

	if streq(f.alg, 'map-adp')
		f.out = [f.dir 't0'];
		t = sprintf('op stack %s float t,out,??.fld', f.out)
		os_run(t)
	end
	xasp = double(fld_read(f.out));
end

if 1
	im plc 2 2
	im(1, xmat, 'xhat matlab'), cbar
	im(2, xasp, 'xhat aspire'), cbar
	im(3, (xasp-xmat)/max(xmat(:)), 'aspire-matlab'), cbar

	t = vcorrcoef(xasp, xmat);
	printm('corr. %g,%g', t, t-1)

	t1 = eql_obj(xmat, G, yi(:), ci(:), ri(:), R, ig.mask);
	t2 = eql_obj(xasp, G, yi(:), ci(:), ri(:), R, ig.mask);

	if im
		im subplot 4
		plot(0:f.niter, t1-t1(1), '-o', 0:f.niter, t2-t1(1), '-x')
		legend('mat', 'asp', 4)
		title 'Objective function'
	end
	% todo: compare to lbfgs!
end
