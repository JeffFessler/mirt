% tql_os_sps_test.m
% compare aspire and matlab versions of T-QPL-OS-SPS
% Copyright Apr 2000	Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'), disp 'setup tql_osps_test'
	if has_aspire && has_mex_jf
		f.dir	= test_dir;
		f.wtf	= [f.dir 't,g.wtf'];
		f.wtr	= strrep(f.wtf, 'wtf', 'wtr');
		f.yi	= [f.dir 'yi.fld'];
		f.bi	= [f.dir 'bi.fld'];
		f.ri	= [f.dir 'ri.fld'];
		f.mask	= [f.dir 'mask.fld'];
	end
	tr_test_setup
prompt
end

if ~isvar('R')
	f.l2b = 16;
	f.nbrs = 4;
	R = Robject(ig.mask, 'edge_type', 'leak', 'offsets', [1 ig.nx], ...
		'type_denom', 'matlab', 'beta', 2^f.l2b);
prompt
end

if ~isvar('Gb'), disp 'make Gb'
	f.nblock = 5;
	Gb = Gblock(G, f.nblock);
prompt
end

%
% matlab iterations
%
if ~isvar('xmat'), disp 'matlab T-QPL-OS-SPS'
	f.niter = 9;
	f.curv = 'pc';
%	f.curv = 'oc';
	f.pixmax = 7;

	xinit = xfbp;
	xmat = tpl_os_sps(xinit(ig.mask), Gb, yi, bi, ri, R, ...
		f.niter, f.pixmax, f.curv);
	xmat = ig.embed(xmat);

	im clf, im(xmat, 'Matlab T-QPL-OS-SPS iterations')
prompt
end

if ~has_aspire || ~has_mex_jf, return, end

%
% aspire iterations
%
if ~isvar('xasp'), disp 'aspire T-QPL-OS-SPS'

	f.init	= [f.dir 'init.fld'];
	f.out	= [f.dir 'out.fld'];
	fld_write(f.init, xinit, 'check', 0)
	if exist(f.out, 'file'), delete(f.out), end

	f.saver	= 'stack,1';
	f.alg	= sprintf('ospsc,%s,%d,%d,1,0', f.curv, f.nblock, sg.na);
	f.penal	= sprintf('%g,quad,%d,-', f.l2b, f.nbrs/4);
	f.method = sprintf('@%d@%s@%s', f.niter-1, f.alg, f.penal);
	f.fitype = ['2z@' f.wtr '@-'];
	f.com = sprintf(['i -chat 0 trpl3 %s %s  %s %s 1 %s 1 %s -' ...
			' %s %s 0 1 %g 0 -'], ...
		f.out, f.init, f.yi, f.bi, f.ri, f.fitype, ...
		f.method, f.saver, f.pixmax);

	os_run(f.com)

	xasp = double(fld_read(f.out));
end

if 1
	t = vcorrcoef(xasp, xmat);
	printf('corr. %g,%g', t, t-1)

	im clf, im(221, xmat, 'xhat matlab'), cbar
	im(222, xasp, 'xhat aspire'), cbar
	im(223, (xasp-xmat)/max(xmat(:)), 'aspire-matlab'), cbar

	t1 = tpl_obj(xmat, G, yi(:), bi(:), ri(:), R, ig.mask);
	t2 = tpl_obj(xasp, G, yi(:), bi(:), ri(:), R, ig.mask);

	if im
		subplot(224)
		plot(0:f.niter-1, t1-t1(1), '-o', 0:f.niter-1, t2-t1(1), '-x')
		xlabel iteration, ylabel '\Phi change', legend('mat', 'asp', 4)
		title(sprintf('T-QPL-OS-SPS, Nsubset=%d', f.nblock))
	end
end
