% tml_os_sps_test.m
% compare aspire and matlab versions of T-ML-OS-SPS
% Copyright Apr 2000, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('xfbp'), printm 'setup tml_os_sps_test'
	if has_aspire && has_mex_jf
		f.dir	= test_dir;
		f.dsc	= [f.dir 't,g.dsc'];
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

if ~isvar('Gb'), printm 'Gb'
	f.nblock = 5;
	Gb = Gblock(G, f.nblock);
end

%
% matlab iterations
%
if ~isvar('xmat'), printm 'matlab T-ML-OS-SPS'
	f.niter = 9;
	f.curv = 'pc';
%	f.curv = 'oc';

	xinit = ig.ones * mean(xfbp(ig.mask));
	xinit = xtrue;
	xmat = tpl_os_sps(xinit(ig.mask), Gb, yi, bi, ri, [], f.niter, inf, f.curv);
	xmat = ig.embed(xmat);

	im clf, im(xmat, 'Matlab T-ML-OS-SPS iterations'), cbar
prompt
end

if ~has_aspire || ~has_mex_jf, return, end

%
% aspire iterations
%
if ~isvar('xasp'), printm 'aspire T-ML-OS-SPS'

	f.init	= [f.dir 'init.fld'];
	f.out	= [f.dir 'out.fld'];
	fld_write(f.init, xinit, 'check', 0)
	if exist(f.out, 'file'), delete(f.out), end

	f.saver	= 'stack,1';
	f.alg	= sprintf('ospsc,%s,%d,%d,1,0', f.curv, f.nblock, sg.na);
	f.penal	= '-';		% unpenalized
	f.method = sprintf('@%d@%s@%s', f.niter-1, f.alg, f.penal);
	f.fitype = ['2z@' f.wtr '@-'];
%	f.fitype = ['2dsc@' f.dsc '@-'];
	f.com = sprintf(['i -chat 0 trpl3 %s %s  %s %s 1 %s 1 %s -' ...
		' %s %s 0 1 1e30 0 -'], ...
		f.out, f.init, f.yi, f.bi, f.ri, f.fitype, ...
		f.method, f.saver);

	os_run(f.com)

	xasp = double(fld_read(f.out));
end

if 1
	t = vcorrcoef(xasp, xmat);
	printf('corr. %g,%g', t, t-1)

	im pl 2 2, im(1, xmat, 'xhat matlab'), cbar
	im(2, xasp, 'xhat aspire'), cbar
	im(3, xasp-xmat, 'aspire-matlab'), cbar

	t1 = tpl_obj(xmat, G, yi(:), bi(:), ri(:), [], ig.mask);
	t2 = tpl_obj(xasp, G, yi(:), bi(:), ri(:), [], ig.mask);

	if im
		subplot(224)
		plot(0:f.niter-1, t1-t1(1), '-o', 0:f.niter-1, t2-t1(1), '-x')
		xlabel 'iteration', ylabel '\Phi change'
		legend('mat', 'asp', 4)
		title(sprintf('T-ML-OS-SPS, Nsubset=%d', f.nblock))
	end
end
