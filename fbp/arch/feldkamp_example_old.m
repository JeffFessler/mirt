% feldkamp_example.m
% example of how to use feldkamp.m for cone-beam CT reconstruction
% Copyright 2004-8-28, Nicole Caparanis, Patty Laskowsky, Taka Masuda,
% and Jeff Fessler, The University of Michigan

if ~isvar('proj'), disp 'proj'
	down = 8;
	nh = 256/down;
	nv = 240/down;
	na = 224/down;
	ds = 1024/nh;
	dt = ds;
	dis_src_det = 949.075;
	dis_iso_det = 408.075;
	dis_src_iso = dis_src_det - dis_iso_det;
%	dis_foc_src = inf; % flat detector panel
	dis_foc_src = 0; % arc detector
	offset_det_h = 0.25; % quarter detector
	offset_det_v = 0.0;
	horiz = ([-(nh-1)/2:(nh-1)/2]' - offset_det_h) * ds;
	verti = ([-(nv-1)/2:(nv-1)/2]' - offset_det_v) * dt;
	printf('rmax=%g', dis_src_iso*sin(atan(max(abs(horiz)) / dis_src_det)))

	ell = [ ...
		[20 0 0  150 150 200 0 0.01]; % 30cm diam "cylinder"
		[80 0 0 50 50 30  0 0.01]; % bone-like inserts
		[0 0 75  40 40 40  0 0.01];
		[0 70 0  30 30 30  0 0.01];
	];

	proj = ellipsoid_proj(ell, horiz, verti, na, ...
		dis_src_iso, dis_iso_det, dis_foc_src);

	nx = 256/down;
	ny = 240/down;
	nz = 200/down;
	dx = 2*down; dy = dx; dz = dx; % 2mm voxels * down-sampling
	xtrue = ellipsoids(nx, ny, nz, ell, dx, dy, dz);

	% cone-beam system geometry, generalized from fan-beam geometry.
	% see ASPIRE users guide under tech. reports on web page for details.
	args = arg_pair('system', NaN, 'nx', nx, 'ny', ny, 'nz', nz, ...
		'nv', nv, 'nh', nh, 'na', na, 'support', 'all', ...
		'orbit', 360, 'orbit_start', 0, ...
		'pixel_size', dx, 'ray_spacing', ds, 'strip_width', 0, ...
		'dis_src_det', dis_src_det, ...
		'dis_iso_det', dis_iso_det, ...
		'dis_foc_src', dis_foc_src, ...
		'offset_source', 0, ...
		'offset_det_h', offset_det_h, ...
		'offset_det_v', offset_det_v);

	pl=330;
	im(pl+1, xtrue, 'x true'), cbar
	im(pl+4, proj, 'true projections'), cbar
	drawnow
prompt
end

% noisy data and estimated line integrals
if ~isvar('li_hat'), disp 'li_hat'
	% noisy data, if blank scan value has been specified.
	if isvar('bi') && isvar('ri')
		yb = bi .* exp(-proj) + ri;
		yi = poisson(yb);
		li_hat = -log((yi-ri) ./ bi);
		li_hat(yi-ri <= 0) = 0; % fix: need something better here...
	else
		li_hat = proj; % noiseless
	end
end

% FDK cone-beam reconstruction
if ~isvar('xfdk'), disp 'fdk'
	mask = true([nx ny nz]);
	xfdk = feldkamp_old(li_hat, 'ramp', mask, args);
end

if 1
	% show results (off-center slices worse than central slice)
	im(pl+2, xfdk, 'FDK recon'), cbar
	im(pl+3, xfdk - xtrue, 'FDK error'), cbar

	subplot(pl+5)
	ix = 1:nx; iy = ceil(ny/2); iz = ceil(nz/2);
	plot(ix, xtrue(ix,iy,iz), '-', ix, xfdk(ix,iy,iz), '--')
	axis([1 nx -0.005 0.025]), legend('true', 'FDK recon', 2)
	title 'middle slice', xlabel 'ix'

	subplot(pl+6)
	iz=1:nz; ix = 1+floor(nx/2); iy = 1+floor(ny/2);
	plot(iz, squeeze(xtrue(ix,iy,iz)), '-', iz, squeeze(xfdk(ix,iy,iz)), '--')
	axis([1 nz -0.005 0.025]), legend('true', 'FDK recon', 2), xlabel 'iz'
	title(sprintf('profile at (ix,iy)=(%g,%g)', ix,iy))

prompt
end
