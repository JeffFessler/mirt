% feldkamp_example.m
% example of how to use feldkamp.m for cone-beam CT reconstruction
% Copyright 2004-8-28, Nicole Caparanis, Patty Laskowsky, Taka Masuda,
% and Jeff Fessler, The University of Michigan
% modified version by Ajay Paidi

if ~isvar('proj')
	down = 1/2;
	nv = 240*down;%240*down;
	nh = 256*down;%256*down;
	na = 224*down;%224*down;
	ds = 1024/nh;
	dt = ds;
	dis_src_det = 949;
	dis_iso_det = 408;
	dis_src_iso = dis_src_det - dis_iso_det;
	dis_foc_src = inf; % flat detector panel
	offset_det_h = 0.25; % quarter detector
	offset_det_v = 0.0;
	horiz = ([-(nh-1)/2:(nh-1)/2]' - offset_det_h) * ds;
	verti = ([-(nv-1)/2:(nv-1)/2]' - offset_det_v) * dt;
	printf('rmax=%g', dis_src_iso*sin(atan(max(abs(horiz)) / dis_src_det)))

    
ell = [ ...
		[0 0 -72 150 200 15 0 5];
		[0 0 -32 150 200 20 0 10];
        [0 0 12 150 200 15 0 10];
		[0 0 52 150 200 20 0 7];
        [0 0 82 150 200 15 0 10];
        ];
    
	
	nx = 240*down;%256*down;
	ny = 256*down;%240*down;
	nz = 112*down;%184*down;
	dx = 2/down; dy = dx; dz = dx;
	x = ellipsoids(nx, ny, nz, ell, dx, dy, dz);

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
    
   
	pl=230;
	im(pl+1, x, 'x'), cbar
	im(pl+2, proj, 'proj'), cbar
	drawnow
prompt
end

% cone-beam reconstruction
mask = true([nx ny nz]);
recon = feldkamp(proj, 'ramp', mask, args);

% show results (off-center slices worse than central slice)
im(pl+4, recon, 'recon'), cbar
im(pl+5, recon - x, 'error'), cbar
ix = 1:nx; iy = ny/2; iz = nz/2;
subplot(236)
plot(ix, x(ix,iy,iz), '-', ix, recon(ix,iy,iz), '--')
axis([1 nx -1 16]), legend('true', 'recon', 2)
title 'MIDDLE SLICE'
iz=12;
subplot(233)
plot(ix, x(ix,iy,iz), '-', ix, recon(ix,iy,iz), '--')
axis([1 nx -1 16]), legend('true', 'recon')
title(sprintf('SLICE %d', iz))
%figure
%subplot(2,1,1)
%ix1 = nx/2;
%iy1 = ny/2;
%iy = 1:ny;
%ix = 1:nx;
%iz = nx/2;
%plot(iy, x(ix1, iy, iz), '-', iy, recon(ix1,iy,iz), '--');
%title('plot of x vs. y on middle slice')
