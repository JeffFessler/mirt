% em3_test_setup.m
%
% create sample image, system matrix, and sinograms for examples
% and testing of 3D Poisson emission maximum likelihood (ML) algorithms
% creates: ig3, ig2, ...
%
% Copyright 2001-07-23, Jeff Fessler, The University of Michigan

%
% true emission image
%
if ~isvar('xtrue')
	if ~isvar('ig3')
		% trick: need dx=dy for 3b to work
		ig3 = image_geom('nx', 32, 'ny', 34, 'nz', 6, 'dx', 1, 'dy', 1);
		ig2 = image_geom('nx', ig3.nx, 'ny', ig3.ny, ...
			'dx', ig3.dx, 'dy', ig3.dy);
	end
	xtrue = [0 0 0	ig3.nx/2-4 ig3.ny/2-4 inf 0 0 1; % cyl
		5 0 0	ig3.nx/15*[1 1 1] 0 0 1;
		4 6 2	ig3.nx/13*[1 1 1] 0 0 1;
		-3 0 0.5 ig3.nx/17*[1 1 1] 0 0 1];
	xtrue = ellipsoid_im(ig3, xtrue, 'oversample', 2);

	im plc 3 3
	im(1, xtrue, 'true emission images'), cbar

	% reconstruction mask (which pixels do we estimate?)
	ig3.mask = ig3.circ(ig3.nx/2-2, ig3.ny/2-2, 0, 0);
	ig2.mask = ig3.mask_or;
	im(2, ig3.mask, 'support mask'), cbar
prompt
end


%
% system matrix G
%
if ~isvar('G'), printm 'make G'
	if 0 % 3d: 3b@
		sg = sino_geom('par', 'nb', 90, 'na', 60, 'dr', ig3.dx, ...
			'strip_width', ig3.dx);

		% caution: 3b system model is plain line integrals
		f.sys_type = sprintf('3b@1,1,1,0,0,0@-2d,%d@-%d,%d,%d', ...
			sg.na, sg.nb, ig3.nz, sg.na);

	else % 2.5d
		sg = sino_geom('par', 'nb', 36, 'na', 30, 'dr', ig3.dx, ...
			'strip_width', ig3.dx);

		G2 = Gtomo2_strip(sg, ig2, 'single', 1);

		if ~isvar('f.wtf')
			f.wtf = [test_dir 'tmp.wtf'];
		end
		if 1 % ugly hack to get matlab sparse matrix into aspire file
			if exist(f.wtf), delete(f.wtf), end
			[ii jj ss] = find(G2.arg.G);
			tmp = find(ig2.mask(:));
			jj = tmp(jj);
			tmp = sparse(ii, jj, ss, sg.nb*sg.na, ig2.nx*ig2.ny);
			wtf_write(f.wtf, tmp, ig2.nx, ig2.ny, sg.nb, sg.na);
		end
		if ~isvar('f.wtr')
			f.wtr = [test_dir 'tmp.wtr'];
		end
		if 1
			if exist(f.wtr, 'file'), delete(f.wtr), end
			os_run(sprintf('wt -chat 0 col2row %s %s', f.wtr, f.wtf))
		end

		f.sys_type = sprintf('2z@%s@-', f.wtr);
	end

	if ~isvar('f.chat'), f.chat = 0; end
	if ~isvar('f.nthread'), f.nthread = 1; end
	G = Gtomo3(f.sys_type, ig3.mask, ig3.nx, ig3.ny, ig3.nz, ...
		'nthread', f.nthread, 'chat', f.chat);
end


%
% noisy measurements
%
if ~isvar('yi'), printm 'make yi'
	rng(0)
	proj = G * xtrue;
	count = 1e6;
	% detector efficiency variations per CTI 931 PET scanner
	ci = exp(0.3 * randn(size(proj)));
	ci = count / sum(ci(:) .* proj(:)) * ci;
	ci = dsingle(ci);
	ytrue = ci .* proj;
	if ~isvar('randpercent')
		randpercent = 10;
	end
	ri = randpercent / 100 * mean(ytrue(:)) * ones(size(ytrue));
	ri = dsingle(ri);
	yi = poisson(ytrue + ri);

	im(4, proj, 'proj: ideal projections'), cbar horiz
	im(5, ytrue, 'ytrue: true projections'), cbar horiz
	im(6, yi, 'yi: noisy projections'), cbar horiz
	clear count randpercent ytrue proj
end


%
% FBP reconstruction
%
if ~isvar('xfbp')
	if streq(f.sys_type, '2z', 2)
		xfbp = em_fbp(sg, ig2, yi, ci, ri);
	else
		p = [1 3 2];
		xfbp = em_fbp(sg, ig2, permute(yi,p), permute(ci,p), permute(ri,p));
	end
	xfbp = max(xfbp, 0);
	im(7, xfbp, 'FBP Reconstruction'), cbar
end
return


%
% save to files if needed
%
if isvar('f.yi')
	fld_write(f.yi, yi, 'check', 0)
end
if isvar('f.ci')
	fld_write(f.ci, ci, 'check', 0)
end
if isvar('f.ri')
	fld_write(f.ri, ri, 'check', 0)
end
if isvar('f.mask')
	fld_write(f.mask, mask, 'check', 0)
end
