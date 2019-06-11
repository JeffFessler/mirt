% tpl_os_sps_test.m
% compare aspire and matlab versions of T-PL-OS-SPS
% Copyright Apr 2000, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'), printm 'setup tpl_os_sps_test'
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

if ~isvar('Gb'), printm 'make Gb'
	f.nblock = 5;
	Gb = Gblock(G, f.nblock);
prompt
end

if ~isvar('R'), printm 'make R'
	f.nbrs = 8;
	f.l2b = 20;
	f.delta = 0.001;
	f.pot = 'huber';
	R = Robject(ig.mask, 'edge_type', 'leak', 'beta', 2^f.l2b, ...
		'potential', f.pot, 'delta', f.delta, 'type_denom', 'aspire');
prompt
	clear xmat xasp
end

%
% matlab iterations
%
if ~isvar('xmat'), printm 'matlab T-PL-OS-SPS'
	f.niter = 11;
	f.pixmax = 7;
	xinit = max(xfbp,0);

	curvs = {'oc', 'pc'};
	for ic = 1:length(curvs)
		tmp = tpl_os_sps(xinit(ig.mask), Gb, yi, bi, ri, R, ...
			1+f.niter, f.pixmax, curvs{ic});
		xmat{ic} = ig.embed(tmp);
	end
	im clf, im(xmat{1}, 'Matlab T-PL-OSPS iterations')
prompt
end

%
% iot iterations
%
if ~isvar('xiot'), printm 'matlab T-PL-IOT'
	xiot = pl_iot(xinit(ig.mask), Gb, {yi, bi, ri}, R, ...
			'dercurv', 'trl', ...
			'niter', f.niter, ...
			'os', 3, ...
			'riter', 2, ... % change to 1 to match OS
			'pixmax', f.pixmax, ...
			'curvtype', curvs{2});
	xiot = ig.embed(xiot);
	im clf, im(xiot, 'Matlab T-PL-IOT iterations')
prompt
end


if ~has_aspire || ~has_mex_jf, return, end

%
% aspire iterations
%
if ~isvar('xasp'), printm 'aspire T-QPL-OSPS'

	f.init	= [f.dir 'init.fld'];
	f.out	= [f.dir 'out.fld'];
	fld_write(f.init, xinit, 'check', 0)

	if f.nbrs == 8
		f.hood_aspire = 6;
	elseif f.nbrs == 4
		f.hood_aspire = 5;
	else
		error hood
	end

	f.fitype = ['2z@' f.wtr '@-'];
	f.saver	= 'stack,1';
	f.penal	= sprintf('3d,%g,-100,%s,%d,-,%g,ih,1', ...
		f.l2b, f.pot, f.hood_aspire, f.delta);
%	f.penal	= sprintf('%g,%s,%d,-,%g,ih,1', f.l2b, f.pot, f.nbrs/4, f.delta);

	for ic = 1:length(curvs)
		f.alg = sprintf('ospsc,%s,%d,%d,1.,0.', curvs{ic}, f.nblock, sg.na);
		f.method = sprintf('@%d@%s@%s', f.niter, f.alg, f.penal);

		if exist(f.out, 'file'), delete(f.out), end
		f.com = sprintf(['i -chat 0 trpl3 %s %s  %s %s 1 %s 1 %s -' ...
			' %s %s 0 1 %g 0 -'], ...
			f.out, f.init, f.yi, f.bi, f.ri, f.fitype, ...
			f.method, f.saver, f.pixmax);

		os_run(f.com)

		xasp{ic} = double(fld_read(f.out));
	end
end

if 1
	for ic = 1:length(curvs)
		t = vcorrcoef(xasp{ic}, xmat{ic});
		printf('corr. %g,%g', t, t-1)
	end

	im clf, im(221, xmat{1}, 'xhat matlab'), cbar
	im(222, xasp{1}, 'xhat aspire'), cbar
	im(223, (xasp{1}-xmat{1})/max(col(xmat{1})), '(aspire-matlab)/max(matlab)'), cbar

	t1 = tpl_obj(xmat{2}, G, yi(:), bi(:), ri(:), R, ig.mask);
	t2 = tpl_obj(xasp{2}, G, yi(:), bi(:), ri(:), R, ig.mask);
	t0 = tpl_obj(xiot, G, yi(:), bi(:), ri(:), R, ig.mask);

	if im
		subplot(224)
		plot(	0:f.niter, t1-t1(1), '-o', ...
			0:f.niter, t2-t1(1), '-x', ...
			0:f.niter, t0-t1(1), '-+')
		xlabel iteration, ylabel '\Phi change'
		legend('mat', 'asp', 'iot', 4)
		title(sprintf('T-PL-OSPS, Nsubset=%d', f.nblock))
	end
end

if 0 % show final images
	im clf, im([xfbp, xasp{2}(:,:,end); xmat{2}(:,:,end), xiot(:,:,end)]), cbar
end
