% helix_example.m
% Example of how to use rebin_helix.m for helical cone-beam CT reconstruction.
% Copyright 2010-7-21, Gregory Handy and Jeff Fessler, University of Michigan

if 0 % test a particular geometry (case 2795)
	f.down = 4; % down sample a lot to save time
	cg = ct_geom('ge2', 'nt', 32, 'na', 3625, ...
...%		'dfs', inf, ...
		'dfs', 0, ...
		'down', f.down, ...
		'pitch', .53125, 'source_z0', -50.5494160, ...
		'orbit', 3625/984*360, ... % 1326.21951219512195,
		'orbit_start', 109.12139);

	ig = image_geom('nx', 320, 'ny', 320, 'nz', 3*64, ... % 61
		'down', f.down, ...
		'dx', 2.191162, 'dz', 0.625, 'offset_z', 56.0);
%	cg.plot(ig), return
end


if ~isvar('cg'), printm 'cg: cone-beam CT geometry'
	f.down = 4; % down sample a lot to save time
	f.nturn = 12;
	cg = ct_geom('ge2', 'pitch', 0.5, 'source_z0', -100, ...
		'na', 984*f.nturn, 'orbit', 360*f.nturn, 'orbit_start', 17, ...
		'down', f.down);
end

if ~isvar('ig'), printm 'ig: image geometry'
	ig = image_geom('nx', 256, 'ny', 256, 'nz', 160, ...
			'dx', 2, 'dz', 0.625, 'down', f.down);
	im clf, cg.plot(ig);
prompt
end


if ~isvar('ell'), printm 'ell: ellipsoid object'
	z0 = ig.offset_z * ig.dz;
	ell = [ ...
		[0 0 -z0	[0.3 0.1]*ig.fov 0.6*ig.zfov	0 0 1000];
		[80 10 -z0	30 30 10	0 0 1000];
		[-10 -40 75	40 40 40	0 0 1000];
		[-10 80 -20	30 30 30	0 0 1000];
		[0 0 -36	20 20 5	0 0 1000];
		[0 0 -12	20 20 5	0 0 1000];
		[0 0 12	20 20 5	0 0 1000];
		[0 0 36	20 20 5	0 0 1000];
	];
	clear z0
end


if ~isvar('xtrue'), printm 'xtrue: true image volume'
	xtrue = ellipsoid_im(ig, ell, 'oversample', 2);

	clim = [0 2000];
	if im
		figure(1), im plc 2 3
%		im(1, ig.x, ig.y, xtrue, clim), cbar
		im(1, 'mid3', xtrue, clim), cbar
		titlef('x true, z=%g to %g', ig.z(1), ig.z(end))
	end
prompt
end


if ~isvar('proj'), printm 'proj: analytical ellipsoid projection views'
	proj = ellipsoid_proj(cg, ell);

	im(4, 'row', 1, permute(proj, [1 3 2]), 'true helix sinograms'), cbar
prompt
end


if ~isvar('sino'), printm 'sino: rebin helix to fan-beam'
	f.short = 1;
	f.itype = 'linear';
%	f.itype = 'nearest';
	[sino, orbits, used] = rebin_helix(cg, ig, proj, ...
		'type', f.itype, 'short', f.short, 'collapse', 0);

	if exist('helix_rebin_mex') == 3 % matlab vs mex
		[sino0, orbits0, used] = rebin_helix(cg, ig, proj, 'use_mex', 0, ...
			'type', f.itype, 'short', f.short, 'collapse', 0);

		equivs(orbits, orbits0)
		equivs(sino0, sino)
		clear sino0 orbits0
	end

	im(5, sino, 'SSRB sinos'), cbar, axis normal
	im(6, used), cbar, axis normal
prompt
end


% helical cone-beam reconstruction based on SSRB
if ~isvar('xssrb'), printm 'fbp from ssrb'
	[xssrb, tmp] = fbp_helix_stack(cg, ig, sino, orbits, ...
		'short', f.short, 'window', 'hanning,1.0');
	xssrb = xssrb .* (ig.circ(cg.rmax) > 0); % mask

%	im(2, xssrb, 'SSRB recon', clim), cbar
	im(2, 'mid3', xssrb, 'SSRB recon', clim), cbar
	im(6, tmp, 'after weighting'), cbar, axis normal
	im(3, xssrb - xtrue, 'error', [-500 500]), cbar
prompt
end


if im % profile, for checking HU accuracy
	im subplot 6
	iz = 1:ig.nz;
	ix = ig.nx/2+1; iy = ig.ny/2+1;
	plot(	ig.z, squeeze(xtrue(ix,iy,iz)), 'o-', ...
		ig.z, squeeze(xssrb(ix,iy,iz)), '.--')
	axis([minmax(ig.z)' -200 2100])
	titlef('profile at (ix,iy)=(%g,%g)', ix,iy), xlabel z
end

if 0 % toggle
	fun = @(x) x;
	fun = @(x) jf_mip3(x, 'type', 'mid');
	figure(2)
	im_toggle(fun(xssrb .* ig.circ), fun(xtrue))
return
end


return % below here is comparisons with GH original method

if ~isvar('xssrb_gh'), printm 'fbp from ssrb gh'
	[xssrb_gh, sino_gh] = fbp_helix_gh(cg, ig, proj, 'short', f.short);
	im(3, xssrb_gh, 'SSRB GH recon', clim), cbar
	im(4, sino, 'SSRB GH sinos'), cbar, axis normal
prompt
end

if 0
	figure(2)
%	im_toggle(fun(xssrb .* ig.circ), fun(xtrue), fun(xssrb_gh .* ig.circ), clim)
	im_toggle(fun(xssrb .* ig.circ), fun(xssrb_gh .* ig.circ), clim)
end

if im % profile
	iz = 21;
	iz = 1:ig.nz;
	im(3, xssrb_gh(:,:,iz) - xtrue(:,:,iz), 'GH Error'), cbar

	im subplot 6
	ix = 1:ig.nx; iy = ceil(ig.ny/2); iz = ceil(ig.nz/2);
	iz = round(25/40 * ig.nz);
	plot(ix, xtrue(ix,iy,iz), '-', ix, xssrb_gh(ix,iy,iz), '--')
	axis([1 ig.nx -200 2100])
	titlef('slice %d', iz), xlabel 'ix'
end

% xfdk = easyhelix(cg, ig, li_hat, 'use_mex', has_mex_jf, 'w1cyl', 0);
% xfdk = myFeldkamp(cg, ig, li_hat, 'use_mex', 0);
% im(4, xfdk(:,:,25), 'FDK recon Other'), cbar
