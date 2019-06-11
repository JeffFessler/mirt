%| pet_2z_example.m
%|
%| Illustrate image reconstruction for multislice 2D PET, aka 2.5D PET,
%| including both intra and inter-slice regularization.
%|
%| Copyright 2005-6-13, Jeff Fessler, University of Michigan

% true emission map
if ~isvar('xtrue'), printm 'xtrue'
	ig = image_geom('nx', 64, 'ny', 56, 'nz', 30, 'dx', 7.5);

	xtrue = ellipsoid_im(ig, ...
		[0 0 0	200 150	400	0 0 1; % body
		100 0 0	20 20 20	0 0 3-2; % lung spot
		0 50 0	30 30 30	0 0 3-2; % heart spot
		], 'oversample', 2);

	im clf, im(ig.x, ig.y, xtrue, 'xtrue')
prompt
end

% (nonuniform) attenuation map
if ~isvar('mumap'), printm 'mumap'
	% attenuation coefficients [1/mm]
	f.mu_water = 0.01;
	f.mu_lung = 0.004;

	mumap = ellipsoid_im(ig, ...
		[0 0 0, 200 150 400	0 0 f.mu_water; % body
		100 0 0, 50 70 90	0 0 f.mu_lung-f.mu_water; % lung
		-100 0 0, 50 70 90	0 0 f.mu_lung-f.mu_water; % lung
		], 'oversample', 2);

	im clf, im(mumap, 'mumap')
prompt

	% support mask (cylindrical)
	ig.mask = ig.circ(ig.dx * (ig.nx/2-2), ig.dy * (ig.ny/2-4)) > 0;
	im clf, im(ig.mask + xtrue + 100 * mumap)
prompt
end

if ~has_mex_jf, printm 'no wtfmex so ending', return, end

% system model
if ~isvar('G'), printm 'G'
	f.dir = test_dir;
	f.dsc = [test_dir 't.dsc'];
	f.wtr = strrep(f.dsc, 'dsc', 'wtr');
	f.mask = [test_dir 'mask.fld'];
	fld_write(f.mask, ig.mask)

	sg = sino_geom('par', 'nb', ig.nx, 'na', 60 * ig.nx / 64, ...
		'dr', ig.dx, 'strip_width', 2*ig.dx);

	if 1 % new way
		tmp = Gtomo2_wtmex(sg, ig, 'mask', ig.mask_or);
		[tmp dum dum dum dum is_transpose] = ...
			wtfmex('asp:mat', tmp.arg.buff, int32(0));
		if is_transpose
			tmp = tmp'; % because row grouped
		end
%		tmp = tmp(:,ig.mask_or);
		delete(f.wtr)
		wtf_write(f.wtr, tmp, ig.nx, ig.ny, sg.nb, sg.na, 'row_grouped', 1)
	else % old way
		t = 'wt -chat 0 dsc -support "file %s" 2 nx %d ny %d nb %d na %d pixel_size %g scale 0 strip_width %g';
		t = sprintf(t, f.mask, ig.nx, ig.ny, sg.nb, sg.na, ig.dx, sg.strip_width)
		os_run([t ' >! ' f.dsc])
		os_run(sprintf('echo y | wt gen %s row', f.dsc)) % row grouped for OS
	end

	f.sys_type = sprintf('2z@%s@-', f.wtr);

	G = Gtomo3(f.sys_type, ig.mask, ig.nx, ig.ny, ig.nz, ...
		'chat', 0, 'view2d', 1, 'nthread', jf('ncore'));

	if 0
		cpu etic
		G * ig.ones;
		cpu etoc
	return
	end
prompt
end

% simulate noisy projection data
if ~isvar('yi'), printm 'yi'
	ytrue = G * xtrue;
	f.counts = 5e5 * ig.nz;
	% attenuation; uniform detector efficiencies for simplicity
	li = G * mumap;
%	minmax(li)
	ci = exp(-li);
	ci = ci * f.counts / sum(col(ci .* ytrue));
	ytrue = ci .* ytrue;
	printf('counts per slice = %g', sum(ytrue(:)) / ig.nz)
	f.scatter_percent = 10;
	% uniform "scatter background" for simplicity
	ri = ones(size(ytrue)) * f.scatter_percent / 100 * mean(ytrue(:));
	yi = poisson(ytrue + ri);
	im(yi, sprintf('yi : %d sinograms', ig.nz))
prompt
end

% FBP recon for comparison
if ~isvar('xfbp'), printm 'fbp'
	sino = permute((yi - ri) ./ ci, [1 3 2]);

	tmp = fbp2(sg, ig, 'window', 'hann');
	xfbp = fbp2(sino, tmp);%, 'window', 'hann');

	im(xfbp, 'fbp'), cbar
	clear sino tmp
prompt
end

% prepare for iterative recon
if ~isvar('Gb'), printm 'Gb'
	f.nblock = sg.na / 4; % lots of subsets!
	Gb = Gblock(G, f.nblock);
end

% check that Gtomo3 works with subsets for 2z
if 0
	y0 = G * xtrue;
	y1 = Gb{1} * xtrue;
	nblock = block_ob(Gb, 'n');
	ia = 1:nblock:sg.na;
	minmax(y1 - y0(:,:,ia))
	im(y1)

	y0 = zeros(size(y0));
	y0(:,:,ia) = y1;
	x0 = G' * y0;
	im(x0)

	x1 = Gb{1}' * y1;
	im(x1)
	minmax(x1 - x0)
return
end

if ~isvar('os_data'), printm 'os_data'
	os_data = {reshaper(yi, '2d'), reshaper(ci, '2d'), ...
		reshaper(ri, '2d')}; % all the data as 2d arrays
end

%
% ML case
%
if 0

	% run OSEM iteration(s) to init!
	if ~isvar('xinit'), printm 'xinit'
		if 1
			xinit = double6(ig.mask);
		else
			xinit = eml_osem(ig.mask(ig.mask), Gb, os_data{:}, ...
					9+1, [], 'classic');
			xinit = ig.embed(xinit);
			xinit = xinit(:,:,:,end);
		end
		im(xinit, 'x init'), cbar
		clear xmlem xosem xinc1 xinc3 like
	prompt
	end

	% ordinary EM
	if ~isvar('xmlem'), printm 'xmlem'
		f.niter = 30+1;
		xmlem = eml_em(xinit(ig.mask), Gb, yi(:), ci(:), ri(:), [], ...
				f.niter);
		xmlem = ig.embed(xmlem);
		im(xmlem(:,:,:,end), 'x ML-EM'), cbar
	prompt
	end

	% OS EM
	if ~isvar('xosem'), printm 'xosem'
		xosem = eml_osem(xinit(ig.mask), Gb, os_data{:}, ...
				f.niter, [], 'classic');
		xosem = ig.embed(xosem);
		im(xosem(:,:,:,end), 'x OS-EM'), cbar
	prompt
	end

	% incremental EM
	if ~isvar('xinc1'), printm 'xinc1'
		f.os1 = 5;
		xinc1 = eml_inc_em(xinit(ig.mask), Gb, os_data{:}, ...
				'niter', f.niter, 'os', f.os1);
		xinc1 = ig.embed(xinc1);
		im(xinc1(:,:,:,end), 'x ML Inc EM'), cbar
	prompt
	end

	if ~isvar('xinc3'), printm 'xinc3'
		xinc3 = eml_inc_em(xinit(ig.mask), Gb, os_data{:}, ...
				'niter', f.niter, 'os', f.os1, 'hds', 3);
		xinc3 = ig.embed(xinc3);
		im(xinc3(:,:,:,end), 'x ML Inc EM 3'), cbar
	prompt
	end

	if ~isvar('like')
		like_fun = @(x) eql_obj(x, G, yi(:), ci(:), ri(:), [], ig.mask);
		like.mlem = like_fun(xmlem);
		like.osem = like_fun(xosem);
		like.inc1 = like_fun(xinc1);
		like.inc3 = like_fun(xinc3);
	end

	if im
		ii = 0:f.niter-1;
		plot( ...
			ii, like.osem, '-o', ...
			ii, like.inc3, '-+', ...
			ii, like.inc1, '-s', ...
			ii, like.mlem, '-x')
		ir_legend({'OSEM', 'INC EM 3', 'INC EM 1', 'MLEM'})
		title 'Log-likelihood vs Iteration'
%		axisy(6.987e7, 6.993e7)
	end


%% penalized-likelihood / MAP case
else
	if ~isvar('R'), printm 'R'
		f.l2b = 2 - 3*log2(ig.nx/64);
		R = Reg1(ig.mask, 'edge_type', 'tight', ...
			'beta', 2^f.l2b, ...
			'type_denom', 'matlab');
		if 0 % predicted PSF for helping choose beta
			wi = ci(:).^2 ./ max(yi(:),1);
			psf = qpwls_psf(G, R, 1, ig.mask, Gdiag(wi));
			im(psf, 'psf'), cbar
		return
		end, clear wi
	end

	% initialize
	if ~isvar('xinit'), printm 'xinit'
		xinit = max(xfbp, 0.1);
		% xinit = double6(ig.mask);
		if 1 % initialize with some OSDP iteration(s) (!)
			xinit = eql_os_emdp(xinit(ig.mask), Gb, os_data{:}, ...
					R, 'niter', 1, 'isave', 1);
			xinit = ig.embed(xinit);
			im(xinit, 'x init'), cbar
		end
		clear xemdp1 xosdp1 xosdp3 xinc1 xinc3 obj
	prompt
	end

	% OSDP1 (non convergent)
	if ~isvar('xosdp1'), printm 'xosdp1'
		f.niter = 30;
		xosdp1 = eql_os_emdp(xinit(ig.mask), Gb, os_data{:}, ...
				R, 'niter', f.niter);
		xosdp1 = ig.embed(xosdp1);
		im(xosdp1(:,:,:,end), 'x OSDP1'), cbar
	prompt
	end

	% EMDP (slow, convergent)
	if ~isvar('xemdp1'), printm 'xemdp1'
		xemdp1 = eql_os_emdp(xinit(ig.mask), Gblock(G,1), ...
				yi(:), ci(:), ri(:), ...
				R, 'niter', f.niter);
		xemdp1 = ig.embed(xemdp1);
		im(xemdp1(:,:,:,end), 'x EMDP1'), cbar
	prompt
	end

	% OSDP3 (non convergent)
	if 0 && ~isvar('xosdp3'), printm 'xosdp3'
		xosdp3 = eql_os_emdp(xinit(ig.mask), Gb, os_data{:}, ...
				R, 'niter', f.niter, 'hds', 3);
		xosdp3 = ig.embed(xosdp3);
		im(xosdp3(:,:,:,end), 'x OSDP3'), cbar
	prompt
	end

	% INC1 (convergent)
	if ~isvar('xinc1'), printm 'xinc1'
		f.os1 = 2;
		xinc1 = epl_inc(xinit(ig.mask), Gb, os_data{:}, ...
			R, 'os', f.os1, 'niter', f.niter, 'hds', 1);
		xinc1 = ig.embed(xinc1);
		im(xinc1(:,:,:,end), 'x inc1'), cbar
	prompt
	end

	% INC3 (convergent)
	if ~isvar('xinc3'), printm 'xinc3'
		xinc3 = epl_inc(xinit(ig.mask), Gb, os_data{:}, R, ...
			'os', f.os1, 'niter', f.niter, 'hds', 3);
		xinc3 = ig.embed(xinc3);
		im(xinc3(:,:,:,end), 'x inc3'), cbar
	prompt
	end

	% test new version
	if 1 && ~isvar('xinc4'), printm 'xinc4'
		xinc4 = epl_inc2(xinit(ig.mask), Gb, yi, ci, ri, R, ...
			'os', f.os1, 'niter', f.niter, 'hds', 3);
		xinc4 = ig.embed(xinc4);
		im(xinc4(:,:,:,end), 'x inc4'), cbar
		minmax(xinc3-xinc4, 'old vs new epl_inc')
	prompt
	end

	if ~isvar('obj')
		obj_fun = @(x) eql_obj(x, G, yi(:), ci(:), ri(:), R.C, ig.mask);
		obj.emdp1 = obj_fun(xemdp1);
		obj.osdp1 = obj_fun(xosdp1);
%		obj.osdp3 = obj_fun(xosdp3);
		obj.inc3 = obj_fun(xinc3);
		obj.inc1 = obj_fun(xinc1);
	end

	if im
		ii = 0:f.niter;
		o0 = obj.emdp1(1);
		plot( ...
			ii(1), -o0+obj.inc3(1), 'b-v', ...
			ii(1), -o0+obj.inc1(1), 'g-x', ...
...%			ii(1), -o0+obj.osdp3(1), 'y-+', ...
			ii(1), -o0+obj.osdp1(1), 'r-o', ...
			ii(1), -o0+obj.emdp1(1), 'c-s', ...
			ii(1:4:end), -o0+obj.inc3(1:4:end), 'bv', ...
			ii(2:4:end), -o0+obj.inc1(2:4:end), 'gx', ...
...%			ii(1:4:end), -o0+obj.osdp3(1:4:end), 'y+', ...
			ii(3:4:end), -o0+obj.osdp1(3:4:end), 'ro', ...
			ii(4:4:end), -o0+obj.emdp1(4:4:end), 'cs', ...
			ii, -o0+obj.inc3, 'b-', ...
			ii, -o0+obj.inc1, 'g-', ...
...%			ii, -o0+obj.osdp3, 'y-', ...
			ii, -o0+obj.osdp1, 'r-', ...
			ii, -o0+obj.emdp1, 'c-')
		ir_legend(	...
			{'C-OS-3', 'INC EM', ...
...%			'OSDP3', ...
			'OSDP', ...
			'EMDP'})
		ylabel 'Penalized-Likelihood Objective'
		xlabel 'Iteration'
		title(mfilename, 'interpreter', 'none')
	end

	printm 'max diffs:'
	disp([max(col(xinc3(:,:,:,end) - xemdp1(:,:,:,end))) ...
	max(col(xinc3(:,:,:,end) - xinc1(:,:,:,end))) ...
	max(col(xinc1(:,:,:,end) - xemdp1(:,:,:,end))) ...
	max(col(xinc3(:,:,:,end) - xosdp1(:,:,:,end))) ...
	max(col(xinc3(:,:,:,end) - xtrue))])
end
