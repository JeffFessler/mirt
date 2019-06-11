%| spect_3d_example.m
%|
%| Illustrate 3D SPECT image reconstruction with compensation
%| for depth-dependent detector response and nonuniform attenuation.
%|
%| For simplicity, there is no "model mismatch" in this example:
%| we use a discrete object and we simulate the projection views
%| with exactly the same system model as is used for reconstruction.
%|
%| The motivation for this example was to illustrate how incremental EM methods
%| (which are convergent) avoid the stagnation of conventional OS methods.  
%| This is also a useful starting point for experimenting with various
%| system models for SPECT reconstruction.
%|
%| Copyright 2005-1-31, Jeff Fessler, University of Michigan

% true emission map
if ~isvar('xtrue'), printm 'xtrue'
	ig = image_geom('nx', 64, 'nz', 30, 'dx', 7.5);

	tmp = [100 0 0	20 20 20	0 0 1]; % lung spot
	voi{1} = ellipsoid_im(ig, tmp, 'oversample', 2);
	tmp = [0 50 0	30 30 30	0 0 1]; % heart spot
	voi{2} = ellipsoid_im(ig, tmp, 'oversample', 2);
	xtrue =	[0 0 0	200 150 400	0 0 1]; % body
	xtrue = ellipsoid_im(ig, xtrue, 'oversample', 2) + voi{1} + voi{2};

	% support mask (cylindrical field of view)
	ig.mask = ig.circ(ig.dx * (ig.nx/2-2), ig.dy * (ig.ny/2-1)) > 0;

	im clf, im('colorneg', xtrue - eps*ig.mask, 'xtrue')
prompt
end

% (nonuniform) attenuation map
if ~isvar('mumap'), printm 'mu map'
	% attenuation coefficients [1/mm]
	f.mu_water = 0.01;
	f.mu_lung = 0.004;
	mumap = ellipsoid_im(ig, ...
		[0 0 0	200 150 400	0 0 f.mu_water; % body
		100 0 0	50 70 90	0 0 f.mu_lung-f.mu_water; % lung
		-100 0 0 50 70 90	0 0 f.mu_lung-f.mu_water; % lung
		], 'oversample', 2);
%	im clf, im('colorneg', mumap - eps*ig.mask, 'mu map')
	im clf, im( mumap, 'mu map'), cbar
prompt
end

% system model
if ~isvar('G'), printm 'G'
	f.dir = test_dir;
	f.file_mumap = [f.dir 'mumap.fld'];
	fld_write(f.file_mumap, mumap); % save mu map to file

	sg = sino_geom('par', 'nb', ig.nx, 'na', 60 * ig.nx / 64, ...
		'orbit', 360, 'dr', ig.dx);

	% make up a simple depth-dependent detector response model
	f.fwhm_collimator = 4;
	f.fwhm_iso = 15; % linearly depth-dependent gaussian blur
	f.blur = sprintf('gauss,%g,%g', f.fwhm_collimator, f.fwhm_iso);

	% note: if you have your own system PSF in file, you can change
	% the following lines (see ASPIRE 3D documentation) to read in your
	% own set of depth-dependent PSF arrays instead of using gaussian
	f.sys_type = sprintf(...
		'3s@%g,%g,%g,%g,%g,1,%s@%s@-@-%d,%d,%d', ...
		ig.dx, abs(ig.dy), abs(ig.dz), sg.orbit, sg.orbit_start, ...
		f.blur, f.file_mumap, sg.nb, ig.nz, sg.na);

	G = Gtomo3(f.sys_type, ig.mask, ig.nx, ig.ny, ig.nz, ...
		'nthread', jf('ncore'), 'chat', 0);

	% for chang attenuation correction we also need a system model
	% that has no detector response in it.

	f.sys_type0 = sprintf(...
		'3s@%g,%g,%g,%g,%g,1,%s@%s@-@-%d,%d,%d', ...
		ig.dx, abs(ig.dy), abs(ig.dz), sg.orbit, sg.orbit_start, ...
		'none', f.file_mumap, sg.nb, ig.nz, sg.na);

	G0 = Gtomo3(f.sys_type0, ig.mask, ig.nx, ig.ny, ig.nz, ...
		'nthread', jf('ncore'), 'chat', 0);
prompt
end

% simulate noisy projection data
if ~isvar('yi'), printm 'yi'
	ytrue = G * xtrue;
	f.counts = 5e5 * ig.nz;
	% uniform "flood map" for simplicity
	ci = ones(size(ytrue)) * f.counts / sum(ytrue(:));
	ytrue = ci .* ytrue;
	printf('counts per slice = %g', sum(ytrue(:)) / ig.nz)
	f.scatter_percent = 10;
	% uniform "scatter background" for simplicity
	ri = ones(size(ytrue)) * f.scatter_percent / 100 * mean(ytrue(:));
	yi = poisson(ytrue + ri);
	im(yi, sprintf('yi : %d projection views', sg.na))
prompt
end


% FBP recon for comparison
if ~isvar('xfbp'), printm 'fbp'
	sino = permute(yi, [1 3 2]); % projections to sinograms

	tmp = fbp2(sg, ig);
	xfbp = fbp2(sino, tmp, 'window', 'hann');
	clear sino tmp

	% chang attenuation correction
	chang = G0' * ones(size(yi));
	chang(chang == 0) = inf;
	chang = sg.na ./ chang;
	chang = chang / 2; % fix: not sure why!?
%	im(chang, 'chang'), cbar
	xfbp = xfbp .* chang;

	im(xfbp, 'fbp'), cbar
prompt
end


% prepare for iterative recon
if ~isvar('Gb'), printm 'Gb'
	f.nblock = sg.na / 4; % lots of subsets!
	Gb = Gblock(G, f.nblock);
end

if ~isvar('os_data'), printm 'os_data'
	os_data = {reshaper(yi, '2d'), reshaper(ci, '2d'), ...
		reshaper(ri, '2d')}; % all the data as 2d arrays
end


% ML case
if 1

	% choose initial image
	if ~isvar('xinit'), printm 'xinit'
		if 1
			xinit = double6(ig.mask); % uniform
		else % osem
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
		f.niter = 20;
		xmlem = eml_em(xinit(ig.mask), Gb, yi(:), ci(:), ri(:), [], ...
				f.niter+1);
		xmlem = ig.embed(xmlem);
		im(xmlem(:,:,:,end), 'x ML-EM'), cbar
	prompt
	end

	% OS EM
	if ~isvar('xosem'), printm 'xosem'
		xosem = eml_osem(xinit(ig.mask), Gb, os_data{:}, ...
				'niter', f.niter);
		xosem = ig.embed(xosem);
		im(xosem(:,:,:,end), 'x OS-EM'), cbar
	prompt
	end

	% incremental EM
	if 0 && ~isvar('xinc1'), printm 'xinc1'
		f.os1 = 5;
		xinc1 = eml_inc_em(xinit(ig.mask), Gb, os_data{:}, ...
				'niter', f.niter, 'os', f.os1);
		xinc1 = ig.embed(xinc1);
		im(xinc1(:,:,:,end), 'x ML Inc EM'), cbar
	prompt
	end

	if 0 && ~isvar('xinc3'), printm 'xinc3'
		xinc3 = eml_inc_em(xinit(ig.mask), Gb, os_data{:}, ...
				'niter', f.niter, 'os', f.os1, 'hds', 3);
		xinc3 = ig.embed(xinc3);
		im(xinc3(:,:,:,end), 'x ML Inc EM 3'), cbar
	prompt
	end

	if ~isvar('like')
		like.mlem = eql_obj(xmlem, G, yi(:), ci(:), ri(:), [], ig.mask);
		like.osem = eql_obj(xosem, G, yi(:), ci(:), ri(:), [], ig.mask);
%		like.inc1 = eql_obj(xinc1, G, yi(:), ci(:), ri(:), [], ig.mask);
%		like.inc3 = eql_obj(xinc3, G, yi(:), ci(:), ri(:), [], ig.mask);
	end

	if im
		ii = 0:f.niter;
		clf, plot( ...
			ii, like.osem, '-o', ...
... %			ii, like.inc3, '-+', ...
... %			ii, like.inc1, '-s', ...
			ii, like.mlem, '-x')
		ir_legend({'OSEM', 'MLEM'})
%		ir_legend({'OSEM', 'INC EM 3', 'INC EM 1', 'MLEM'})
		title 'Log-likelihood vs Iteration'
%		axisy(6.921e7, 6.923e7)
%		axisy(4.5908e7, 4.593e7)
		if 1
			axes('position', [0.2 0.3 0.2 0.2])
			iz = 15;
			ii = f.nblock;
			t = sprintf('EM at %d', ii);
			im(xmlem(:,:,iz,ii+1), t), cbar

			axes('position', [0.5 0.3 0.2 0.2])
			ii = 1;
			t = sprintf('OSEM at %d', ii);
			im(xosem(:,:,iz,ii+1), t), cbar
		end
	prompt
	end

	if im
		clear voi_mlem voi_osem
		for ll=1:2
			voi_true(ll) = sum(col(voi{ll} .* xtrue));
			voi_mlem(:,ll) = ig.maskit(xmlem)' * ig.maskit(voi{ll});
			voi_osem(:,ll) = ig.maskit(xosem)' * ig.maskit(voi{ll});
			voi_mlem(:,ll) = voi_mlem(:,ll) / voi_true(ll) * 100;
			voi_osem(:,ll) = voi_osem(:,ll) / voi_true(ll) * 100;
		end
		ii = 0:f.niter;
		clf, plot(ii, voi_osem, 'c-o', ii, voi_mlem, 'y-x')
		xlabel 'iteration', ylabel '% recovery'
		ir_legend({'OSEM voi1', 'OSEM voi2', 'MLEM voi1', 'MLEM voi2'})
	end


% penalized-likelihood case
else
	if ~isvar('R'), printm 'R'
		f.l2b = 2 - 3*log2(ig.nx/64);
		R = Reg1(ig.mask, 'edge_type', 'tight', ...
			'beta', 2^f.l2b, 'type_denom', 'matlab');
		if 0 % predicted PSF for helping choose beta
			wi = ci(:).^2 ./ max(yi(:),1);
			psf = qpwls_psf(G, R, 1, ig.mask, spdiag(wi));
			im(psf, 'psf'), cbar
		prompt
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
	if 0&~isvar('xosdp3'), printm 'xosdp3'
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
		xinc3 = epl_inc(xinit(ig.mask), Gb, os_data{:}, ...
			R, 'os', f.os1, 'niter', f.niter, 'hds', 3);
		xinc3 = ig.embed(xinc3);
		im(xinc3(:,:,:,end), 'x inc3'), cbar
	prompt
	end

	if ~isvar('obj')
		obj.emdp1 = eql_obj(xemdp1, G, yi(:), ci(:), ri(:), R, ig.mask);
		obj.osdp1 = eql_obj(xosdp1, G, yi(:), ci(:), ri(:), R, ig.mask);
%		obj.osdp3 = eql_obj(xosdp3, G, yi(:), ci(:), ri(:), R, ig.mask);
		obj.inc3 = eql_obj(xinc3, G, yi(:), ci(:), ri(:), R, ig.mask);
		obj.inc1 = eql_obj(xinc1, G, yi(:), ci(:), ri(:), R, ig.mask);
	end

	if 1
		ii = 0:f.niter;
		o0 = obj.emdp1(1);
		plot( ...
			ii(1:6:end), -o0+obj.inc3(1:6:end), 'b-v', ...
			ii(2:6:end), -o0+obj.inc1(2:6:end), 'g-x', ...
...%			ii(1:6:end), -o0+obj.osdp3(1:6:end), '-+', ...
			ii(3:6:end), -o0+obj.osdp1(3:6:end), 'r-o', ...
			ii(4:6:end), -o0+obj.emdp1(4:6:end), 'c-s', ...
			ii, -o0+obj.inc3, '-', ...
			ii, -o0+obj.inc1, '-', ...
...%			ii, -o0+obj.osdp3, '-+', ...
			ii, -o0+obj.osdp1, '-', ...
			ii, -o0+obj.emdp1, '-')
		ir_legend(	...
			{'C-OS-3', 'INC EM', ...
...%			'OSDP3', ...
			'OSDP', ...
			'EMDP'})
		ylabel 'Penalized-Likelihood Objective'
		xlabel 'Iteration'
%		axisy(6.1e5, 6.3e5)
		axisx(0, 190)
	end
	[max(col(xinc3(:,:,:,end) - xemdp1(:,:,:,end))) ...
	max(col(xinc3(:,:,:,end) - xinc1(:,:,:,end))) ...
	max(col(xinc1(:,:,:,end) - xemdp1(:,:,:,end))) ...
	max(col(xinc3(:,:,:,end) - xosdp1(:,:,:,end))) ...
	max(col(xinc3(:,:,:,end) - xtrue))]
end
