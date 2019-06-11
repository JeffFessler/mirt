%| Gtomo3_3s_test.m
%| Test the Gtomo3 object for 3d (SPECT) with many threads vs 1 thread
%|
%| 2017-07-11 updated to test new capability of multiple objects

% f3d_mex('chat', int32(2)) % for debugging

if 0 || ~isvar('x'), printm 'x'
	ig = image_geom('nx', 16, 'ny', 16, 'nz', 10, 'dx', 4, 'dz', 4);
	ig.mask = ig.circ > 0;
	x = single(convn(double(ig.mask), ones(3,3,1)/3^2, 'same') >= 1);
	x(end/4, end/2, 3) = 2;
	im plc 2 3
	im(1, ig.mask)
	im(4, x)
end


if 0 || ~isvar('A1'), printm 'A1'
%	f.chat = 998; % debug
	f.chat = 0;
	f.option = {};

	f3d_mex('free:all') % for testing, want only these two in play

	if 1
		if 1
			f.fwhm_collimator = 1;
			f.fwhm_iso = 2; % depth-dependent gaussian blur
			f.psfs = '-';
			f.blur = sprintf(',gauss,%g,%g', ...
				f.fwhm_collimator, f.fwhm_iso);
		elseif 0
			f.psfs = '-';
			f.blur = ',none';
		else % stress fftw
			f.psfs = '/tmp/t,psfs.fld';
			psfs = make_3s_psfs(ny, 1, 1.2*nx, 0, 2/nx);
			f.blur = ',fft'; % fails due to fftw issues?
			fld_write(f.psfs, psfs)
		end
%		mask = []; % for 3s object!
		f.na = 6;
		f.mumap = '-';
		f.sfilter = 1;
		dx = 4;
		f.sys_type = '3s@%g,%g,%g,360,0,%d%s@%s@%s@-%d,%d,%d';
		f.sys_type = sprintf(f.sys_type, ig.dx, ig.dx, ig.dz, ...
			f.sfilter, f.blur, f.mumap, f.psfs, ig.nx, ig.nz, f.na);

%		3s@[-|sx,sy,sz,orbit,orbit_start,spline_filter[,blur_method]]
%			@mumap.file@filter.file@-nu,nv,nview
	end

	f.option = {f.option{:}, 'chat', f.chat, 'checkmask', 0};

	A1 = Gtomo3(f.sys_type, ig.mask, ig.nx, ig.ny, ig.nz, ...
		f.option{:}, 'nthread', 1);

	tester_tomo2(A1, ig.mask) % put it through paces
end


if 1 || ~isvar('y1'), printm 'y1'
	cpu etic
	y1 = A1 * x;
	t1 = cpu('etoc', sprintf('%d threads', A1.arg.nthread));
	im(2, y1)
end


if 0 || ~isvar('A2'), printm 'A2'
	A2 = Gtomo3(f.sys_type, ig.mask, ig.nx, ig.ny, ig.nz, ...
		f.option{:}, ...
			'nthread', 2+2*jf('ncore'), ...
			'nthread_max', 100);

	tester_tomo2(A2, ig.mask, 'A2', A1) % put it through paces / compare
end


if 0 || ~isvar('y2'), printm 'y2'
	y2 = A2 * x;
	t2 = cpu('etoc', sprintf('%d threads', A2.arg.nthread));

	im(5, y2)
	im(6, y2-y1)

        jf_equal(y1, y2)
end

f3d_mex('show:all')
if 1 % clean up after testing
	clear A1 A2
	f3d_mex('free:all')
end
