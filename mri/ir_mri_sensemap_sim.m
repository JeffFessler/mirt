 function [smap, x, y, z] = ir_mri_sensemap_sim(varargin)
%function [smap, x. y, z] = ir_mri_sensemap_sim(varargin)
%|
%| Simulate 2D or 3D sensitivity maps for sensitivity-encoded MRI
%| based on grivich:00:tmf doi:10.1119/1.19461
%| This code makes maps for multiple coils, but does not model coupling
%| between coils so most likely it is an approximation at best.
%|
%| option
%|	nx, ny, nz		image size (default: [64 64 1])
%|	dx, dy, dz		pixel/voxel dimensions (default: [3 3 3])
%|	ncoil			# of coils total (default: 4)
%|	nring			# of rings of coils (default: 1)
%|	rcoil			coil radius (default: 100mm)
%|	dz_coil			ring spacing in z.  (def: nz*dz/nring)
%|				(3D geometry is a cylinder)
%|	coil_distance		distance of coil center from isocenter for
%|				central ring of coils as a multiple of FOVx,
%|				where FOVx=nx*dx (default: 1.2)
%|	orbit			default: 360
%|	scale			'' (default)
%|				'ssos_center' : make SSoS of center = 1
%|
%| out
%|	smap	[nx ny nz ncoil]	simulated sensitivity maps (complex!)
%|
%| all length parameters must have same units (e.g., mm or cm)
%|
%| Copyright 2005-6-20, Jeff Fessler and Amanda Funai, University of Michigan
%| 2014-08-19 JF more testing, verifying phase is correct
%| 2014-09-09 modified for 3D by Mai Le
%| 2016-05-03 JF fixes 

if nargin < 1, ir_usage, end

if nargin == 1 && streq(varargin{1}, 'test', 4)
	ir_mri_sensemap_sim_test(varargin{1})
	return
end

arg.nx = 64;
arg.ny = [];
arg.nz = 1; % 2D
arg.dx = 3; % pixel size in mm
arg.dy = [];
arg.dz = [];
arg.ncoil = 4; % # of coils
arg.nring = 1;
arg.rcoil = 100; % coil radius in mm
arg.orbit = 360;
arg.orbit_start = 0; % can be [nring] to give each ring an offset [degrees]
arg.dz_coil = [];
arg.coil_distance = 1.2; % multiplies fov/2
arg.scale = '';
arg.chat = nargout == 0;

arg = vararg_pair(arg, varargin);

if isempty(arg.dy), arg.dy = arg.dx; end
if isempty(arg.dz), arg.dz = arg.dx; end
if isempty(arg.ny), arg.ny = arg.nx; end
if isempty(arg.rcoil), arg.rcoil = arg.dx * arg.nx / 2 * 0.50; end
if isempty(arg.dz_coil), arg.dz_coil = arg.dz * arg.nz / arg.nring; end

coils_per_ring = round(arg.ncoil / arg.nring);
if arg.nring * coils_per_ring ~= arg.ncoil
	fail('nring must be divisor of ncoil')
end

[ring_smap, x, y, z] = ir_mri_sensemap_sim_do(...
	arg.nx, arg.ny, arg.nz, ...
	arg.dx, arg.dy, arg.dz, ...
	arg.ncoil, coils_per_ring, arg.rcoil, arg.dz_coil, ...
	arg.orbit, arg.orbit_start, arg.coil_distance, arg.chat);

if arg.nz == 1
	smap = reshape(ring_smap, [arg.nx arg.ny arg.ncoil]);
	scale_center = 1 / sqrt(sum(abs(smap(end/2,end/2,:).^2)));
else
	smap = reshape(ring_smap, [arg.nx arg.ny arg.nz arg.ncoil]);
	scale_center = 1 / sqrt(sum(abs(smap(end/2,end/2,end/2,:).^2)));
end

switch arg.scale
case ''
case 'ssos_center'
	smap = smap * scale_center;
otherwise
	fail('unknown scale method "%s"', arg.scale)
end


% ir_mri_sensemap_sim_do()
function [smap x y z] = ir_mri_sensemap_sim_do(nx, ny, nz, ...
		dx, dy, dz, ncoil, ncoilpr, rcoil, dz_coil, ...
		orbit, orbit_start, coil_distance, chat)

nring = ncoil / ncoilpr;
rlist = rcoil * ones(ncoilpr,nring,'single'); % coil radii

plist = zeros(ncoilpr,nring,3,'single'); % position of coil center [x y z]
nlist = zeros(ncoilpr,nring,3,'single'); % normal vector (inward) from coil center
olist = zeros(ncoilpr,nring,3,'single'); % unit vector orthogonal to normal vector in x-y
ulist = zeros(ncoilpr,nring,3,'single'); % upward vector

if numel(orbit_start) == 1
	orbit_start = repmat(orbit_start, nring);
end

% cylindrical coil configuration, like abdominal coils
alist = deg2rad(orbit) * [0:(ncoilpr-1)] / ncoilpr; % coil angles [radians]
z_ring = ([1:nring]-(nring+1)/2) * dz_coil;
for ir = 1:nring
	for ic = 1:ncoilpr
		phi = alist(ic) + deg2rad(orbit_start(ir));
		Rad = max(nx/2 * dx, ny/2 * dy) * coil_distance;
		plist(ic,ir,:) = [Rad * [cos(phi) sin(phi)] z_ring(ir)];
		nlist(ic,ir,:) = -[cos(phi) sin(phi) 0*z_ring(ir)]; % cylinder
		olist(ic,ir,:) = [-sin(phi) cos(phi) 0];
		ulist(ic,ir,:) = [0 0 1];
	end
end

% object coordinates
x = ([1:nx] - (nx+1)/2) * dx;
y = ([1:ny] - (ny+1)/2) * dy;
z = ([1:nz] - (nz+1)/2) * dz;
[xx,yy,zz] = ndgrid(x,y,z);

smap = zeros(nx, ny, nz, ncoilpr, nring, 'single');
for ir = 1:nring
	for ic=1:ncoilpr
		% rotate coordinates to correspond to coil orientation
		zr =	(xx - plist(ic,ir,1)) .* nlist(ic,ir,1) + ...
			(yy - plist(ic,ir,2)) .* nlist(ic,ir,2) + ...
			(zz - plist(ic,ir,3)) .* nlist(ic,ir,3);
		xr =	xx .* nlist(ic,ir,2) - yy .* nlist(ic,ir,1);
		yr = zz - plist(ic,ir,3); % translate along object z axis

		if 0 % see coordinates
			im plc 1 2
			im(1, x, y, xr), xlabel x, ylabel y
			im(2, x, y, zr)
			keyboard
		end

		% compute sensitivity vectors in coil coordinates
		[sx,sy,sz] = ir_mri_smap1(xr, yr, zr, rlist(ic,ir));

		% coil response depends on tranverse magnetization only?
		% todo: unsure if this should depend on sy and ulist in 3D
		bx = sz * nlist(ic,ir,1) + sx * olist(ic,ir,1);
		by = sz * nlist(ic,ir,2) + sx * olist(ic,ir,2);
	%	bz = sz * nlist(ic,ir,3) + sx * olist(ic,ir,3);
		smap(:,:,:,ic,ir) = bx + 1i * by;

		if 0 && chat && nz == 1 && im % see field components
			im plc 2 2
			im(1, x, y, sx), cbar
			im(2, x, y, sy), cbar
			im(3, x, y, sz), cbar
			im subplot 4
			tmp = sqrt(sx.^2 + sz.^2);
			quiver(x, y, (sx./tmp)', (sz./tmp)', 0), axis square
			prompt
		end

		if 0 % see final field components vs phase
			im subplot 4
			bb = sqrt(bx.^2 + by.^2);
			quiver(x, y, (bx./bb)', (by./bb)', 0), axis square
			im(2, x, y, angle(smap(:,:,ic))), cbar
			keyboard
		end
	end
end
smap = smap * rlist(1) / (2*pi); % trick: scale so maximum is near unity

if chat && im % show smap and array geometry in z=0 plane
	if nz == 1
		ir_mri_sensemap_sim_show2(smap, x, y, dx, dy, nlist, plist, rlist)
	else
		ir_mri_sensemap_sim_show3(smap, x, y, z, dx, dy, dz, ...
		nlist, plist, rlist, olist, ulist, nring, ncoilpr, rcoil)
	end
end


% ir_mri_sensemap_sim_show3()
% shows coil geometry but not the 3D smap
function ir_mri_sensemap_sim_show3(smap, x, y, z, dx, dy, dz, ...
	nlist, plist, rlist, olist, ulist, nring, ncoilpr, rcoil)

pcolor = {'c', 'g', 'r'};
pcolor = @(i) pcolor{1+rem(i,3)};
clf, ir_plot3_cube(x,y,z)
xlabel x, ylabel y, zlabel z
hold on
plot3(plist(:,:,1), plist(:,:,2), plist(:,:,3), 'bo') % coil centers
if 1 % coil normals
	tmp1 = reshape(plist, [], 3);
	tmp2 = reshape(nlist, [], 3);
	quiver3(tmp1(:,1), tmp1(:,2), tmp1(:,3), ...
		tmp2(:,1), tmp2(:,2), tmp2(:,3), 0.2)
end
if 1 % coils
	for ir = 1:nring
	for ic = 1:ncoilpr
		tmp = linspace(0, 2*pi, 50)';
		tmp = cos(tmp) * squeeze(olist(ic,ir,:))' + ...
			sin(tmp) * squeeze(ulist(ic,ir,:))';
		tmp = repmat(squeeze(plist(ic,ir,:))', ...
			[nrow(tmp) 1]) + rcoil * tmp;
	%	plot3(tmp(:,1), tmp(:,2), tmp(:,3), 'g-')
		patch(tmp(:,1), tmp(:,2), tmp(:,3), pcolor(ir), ...
		'edgecolor', 'none', 'facealpha', 0.5)
	end
	end
end
hold off
axis equal
% end


function ir_plot3_cube(x,y,z)
x1 = x(1);
x2 = x(end);
y1 = y(1);
y2 = y(end);
z1 = z(1);
z2 = z(end);
x = [x1 x2 x2 x1 x1 x1 x2 x2 x1 x1];
y = [y1 y1 y2 y2 y1 y1 y1 y2 y2 y1];
z = [z1 z1 z1 z1 z1 z2 z2 z2 z2 z2];
plot3(x,y,z)


% ir_mri_sensemap_sim_show2()
function ir_mri_sensemap_sim_show2(smap, x, y, dx, dy, nlist, plist, rlist)
switch ndims(smap) 
case 3
	[nx ny ncoil] = size(smap);
case 4
	[nx ny nz ncoil] = size(smap);
otherwise
	fail('unknown ndims(smap) = %d', ndims(smap))
end
nshow = min(max(ncoil,2), 4);
im('plc', 3, nshow)
for ii=1:min(ncoil,nshow)
	clim = [0 max(abs(smap(:)))];
	tmp = smap(:,:,ii);
	im(ii, x, y, abs(tmp), clim, 'Magnitude'), cbar
%	xmax = max(max(abs(x)), max(plist(:,1)));
%	ymax = max(max(abs(y)), max(plist(:,2)));
	xmax = max([max(abs(x)) max(abs(y)) max(col(plist(:,[1 2])))]);
%	axis([-xmax xmax -ymax ymax]*1.05)
	axis(xmax * [-1 1 -1 1] * 1.1)
	xtick([-1 0 1] * nx/2 * dx)
	ytick([-1 0 1] * ny/2 * dy)

	hold on
	plot(0,0,'.', plist(:,1), plist(:,2), 'bo')
	xdir = nlist(ii,2);
	ydir = nlist(ii,1);
	r = rlist(ii);
	plot(plist(ii,1)+r*xdir*[-1 1], plist(ii,2)+r*ydir*[1 -1], 'b-')
	hold off

%	ph = ir_unwrap(angle(tmp)); % trick: unwrap phase for pretty display
	ph = angle(tmp); % show raw phase (understandable with hsv colormap)
	im(ii+nshow, 'hsv', x, y, ph, [-pi pi], 'Phase'), cbar
	axis(xmax*[-1 1 -1 1]*1.1)
	xtick([-1 0 1] * nx/2 * dx)
	ytick([-1 0 1] * ny/2 * dy)
end

ssos = sqrt(sum(abs(smap).^2, ndims(smap)));
minmax(ssos)
ssos = ssos / ssos(end/2,end/2);
im(2*nshow+1, x, y, ssos, 'SSoS (normalized)'), cbar
xtick([-1 0 1] * nx/2 * dx)
ytick([-1 0 1] * ny/2 * dy)

if ncoil == 1
	subplot(122)
	bx = real(smap);
	by = imag(smap);
	quiver(x, y, bx', by'), title 'Field pattern in x-y plane'
	axis equal, axis tight
end

% clf, im(angle(smap(:,:,2))), cbar


% ir_mri_smap_r(r, z)
% function for testing near 0
function out = ir_mri_smap_r(r, z)
M = 4 * r ./ ((1 + r).^2 + z.^2); % = k^2, see ellipke
[K E] = ellipke(M);
out = 2 * z ./ r .* ((1 + r).^2 + z.^2).^(-0.5) .* ...
	((1 + r.^2 + z.^2) ./ ((1 - r).^2 + z.^2) .* E - K);


% ir_mri_smap1()
% based on grivich:00:tmf
% for a circular coil in "x-y plane" of radius "a"
% note that coil x-y plane is not same as object x-y plane!
% returns (i,j,k) components of B vector for each (x,y,z) location
function [smap_x smap_y smap_z] = ir_mri_smap1(x, y, z, a)
x = x ./ a; % normalized units
y = y ./ a;
z = z ./ a;
r = sqrt(x.^2 + y.^2);
M = 4 * r ./ ((1 + r).^2 + z.^2); % = k^2, see ellipke
[K E] = ellipke(M);
if ir_is_octave
	K = reshape(K, size(M)); % ellipke shape bug in octave
	E = reshape(E, size(M));
end

% the following is B_z in eqn (18) in grivich:00:tmf
% and same as eqn [10] in wang:00:dop to within constant scale factor
smap_z = 2 * ((1 + r).^2 + z.^2).^(-0.5) .* ...
	(K + (1 - r.^2 - z.^2) ./ ((1 - r).^2 + z.^2) .* E);
smap_z = smap_z / a;

if 0 && any(r(:) == 0) % test code to explore when r is near 0
	r0 = linspace(0,5e-7,101);
	z0 = 0.4;
	t0 = ir_mri_smap_r(r0, z0);
	slope = 3*pi * z0 / ((1+z0^2)^2.5);
	clf, plot(r0, t0, '-', r0, slope * r0, '--'); grid, prompt
end

% the following is B_r in eqn (17) in grivich:00:tmf
smap_r = 2 * z ./ r .* ((1+r).^2 + z.^2).^(-0.5) .* ...
	((1 + r.^2 + z.^2) ./ ((1-r).^2 + z.^2) .* E - K);
bad = abs(r) < 1e-6;
smap_r(bad) = 3 * pi * z(bad) ./ ((1 + z(bad).^2).^2.5) .* r(bad);
smap_r = smap_r / a;

if any(isnan(smap_r(:))) || any(isnan(smap_z(:)))
	keyboard
end

smap_x = smap_r .* div0(x, r);
smap_y = smap_r .* div0(y, r);
%smap_z = smap_r .* div0(z, r);

%phi = atan2(y, x);
%smap_x = smap_r .* cos(phi);
%smap_y = smap_r .* sin(phi);


% ir_mri_sensemap_sim_test0
% see ellipke
function ir_mri_sensemap_sim_test0
m = linspace(0,1,101);
[k e] = ellipke(m);
clf, plot(m, k, '-', m, e, '--'), legend('k', 'e')
yaxis_pi('0 p/2 p 3*p/2')


% ir_mri_sensemap_sim_test1
% test ir_mri_smap1 routine, cf Fig. 4 of grivich:00:tmf
function ir_mri_sensemap_sim_test1
a = 1;
x = linspace(-2,2,99);
y = linspace(-2,2,97);
zlist = [0.001 0.1 0.2 0.5 1.0];
zlist(1) = [];
[xx yy zz] = ndgrid(x, y, zlist);
[smap_x smap_y smap_z] = ir_mri_smap1(xx, yy, zz, a);
smap_b = sqrt(smap_x.^2 + smap_y.^2);
if im
	im('plc', 4, numel(zlist))
	ir_mri_sensemap_sim_test1_show(smap_x, x, y, 0, zlist, 'x')
	ir_mri_sensemap_sim_test1_show(smap_y, x, y, 1, zlist, 'y')
	ir_mri_sensemap_sim_test1_show(smap_z, x, y, 2, zlist, 'z')
	ir_mri_sensemap_sim_test1_show(smap_b, x, y, 3, zlist, 'b')
end


% ir_mri_sensemap_sim_test1_show()
function ir_mri_sensemap_sim_test1_show(map, x, y, offset, zlist, titl)
clim = [-20 20];
for iz = 1:numel(zlist)
	p = offset*numel(zlist) + iz;
	im(p, x, y, map(:,:,iz), titl, clim), cbar
	axis equal
	if iz == 1
		xtick([-2 2]), ytick([-2 2])
	else
		xtick off, ytick off
	end
	if zlist(iz) < 0.5
		blim = [7 12 19];
	else
		blim = [1 3 5];
	end
	hold on
	contour(x, y, abs(map(:,:,iz))', blim, 'b-')
	contour(x, y, abs(map(:,:,iz))', [0 0]+0.001, 'g-')
	hold off
	if streq(titl, 'b')
		xlabelf('z = %g', zlist(iz))
	end
end
drawnow


% ir_mri_sensemap_sim_test2
% 2D test case
function ir_mri_sensemap_sim_test2

[smap x y] = ir_mri_sensemap_sim('chat', 1, 'nx', 32, ...
	'rcoil', [], 'ncoil', 4, 'coil_distance', 1.2);

if 0 % check vs old version
	old = mri_sensemap_sim('chat', 1, 'nx', 32, ...
		'rcoil', [], 'ncoil', 4, 'coil_distance', 1.2);
	equivs(smap, old)
end

if 1 % check rotational symmetry in 4-coil case
	for ic=2:4
		tmp = rot90(smap(:,:,1), (ic-1));
		equivs(abs(tmp), abs(smap(:,:,ic)))
	%	im(8+ic, x, y, abs(tmp)), cbar
		p1 = angle(tmp) + (ic-1) * pi/2; % add pi/2 to rotated
		p2 = angle(smap(:,:,ic));
		tmp = exp(1i * (p2 - p1));
		equivs(tmp, ones(size(tmp))) % trick: equivs mod 2*pi
	end
end


% ir_mri_sensemap_sim_test3
% illustrate 3D sense maps with el 
function ir_mri_sensemap_sim_test3
nring = 3;
ncoil = 4 * nring;
ig = image_geom('nx', 16, 'ny', 14, 'nz', 10, 'fov', 200, 'dz', 20); % 20cm fov
%ig = image_geom('nx', 72, 'ny', 48, 'nz', 12, 'fov', 22, 'zfov', 10); % michelle
%nring = 2; ncoil = 8; % michelle
ig.mask = ig.circ > 0;
smap = ir_mri_sensemap_sim('chat', 1, 'nx', ig.nx, 'ny', ig.ny, 'nz', ig.nz, ...
	'dx', ig.dx, ...
	'dz', ig.dz, ...
	'orbit_start', 1*[0 45 0], ...
...%	'orbit_start', 0*[0 0], ...
	'rcoil', 70, ...
...%	'rcoil', 3, ...
	'nring', nring, 'ncoil', ncoil, 'coil_distance', 1.2);
if im
	prompt
	tmp = smap .* repmat(ig.mask, [1 1 1 ncoil]);
	im clf, im('row', ncoil, reshape(abs(tmp), [ig.nx ig.ny ig.nz*ncoil]))
	prompt
	tmp = permute(tmp, [1 3 2 4]); % [nx nz ny ncoil] z cuts are smooth
	im('row', ncoil, reshape(abs(tmp), [ig.nx ig.nz ig.ny*ncoil]))
	prompt
	im('row', ncoil, abs(tmp(:,:,end/2,:)))
end


% ir_mri_sensemap_sim_test()
function ir_mri_sensemap_sim_test(arg)
switch(arg)
case 'test0'
	ir_mri_sensemap_sim_test0 % ellipk
case 'test1'
	ir_mri_sensemap_sim_test1 % basic test
case 'test2'
	ir_mri_sensemap_sim_test2 % 2D test
case 'test3'
	ir_mri_sensemap_sim_test3 % 3D test
case 'test'
	ir_mri_sensemap_sim_test1
	if im, prompt, end
	ir_mri_sensemap_sim_test2
	if im, prompt, end
	ir_mri_sensemap_sim_test3
otherwise
	fail('bad argument "%s"', arg)
end
