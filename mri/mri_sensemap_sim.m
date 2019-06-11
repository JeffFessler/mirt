 function [smap x y] = mri_sensemap_sim(varargin)
%function [smap x y] = mri_sensemap_sim(varargin)
%|
%| Simulate 2D sensitivity maps for sensitivity-encoded MRI
%| based on grivich:00:tmf doi:10.1119/1.19461
%|
%| option
%|	nx, ny, dx, dy, ncoil, rcoil, coil_distance, orbit (see below)
%|	chat	0|1	show images?
%|
%| out
%|	smap	[nx ny ncoil]	simulated sensitivity maps (complex!)
%|
%| See also ir_mri_sensemap_sim.m for new 3D version!
%|
%| Copyright 2005-6-20, Jeff Fessler and Amanda Funai, University of Michigan
%| 2014-08-19 JF more testing, verifying phase is correct

if nargin == 1 % tests
	if streq(varargin{1}, 'test')
		mri_sensemap_sim_test
	elseif streq(varargin{1}, 'test1')
		mri_sensemap_sim_test1
	end
return
end

arg.nx = 64;
arg.ny = [];
arg.dx = 3; % pixel size in mm
arg.dy = [];
arg.ncoil = 4; % # of coils
arg.rcoil = 100; % coil radius
arg.orbit = 360;
arg.coil_distance = 1.2; % multiplies fov/2
arg.flag_old = false; % old way based on wang:00:dop
arg.chat = nargout == 0;

arg = vararg_pair(arg, varargin);

if isempty(arg.dy), arg.dy = arg.dx; end
if isempty(arg.ny), arg.ny = arg.nx; end
if isempty(arg.rcoil), arg.rcoil = arg.dx * arg.nx / 2 * 0.50; end

[smap x y] = mri_sensemap_sim_do(arg.nx, arg.ny, arg.dx, arg.dy, ...
	arg.ncoil, arg.rcoil, arg.orbit, arg.coil_distance, ...
	arg.flag_old, arg.chat);

if ~nargout, clear, end


% mri_sensemap_sim_do()
function [smap x y] = mri_sensemap_sim_do(nx, ny, dx, dy, ncoil, rcoil, ...
		orbit, coil_distance, flag_old, chat)

rlist = rcoil * ones(ncoil, 1); % coil radii

plist = zeros(ncoil,3); % position of coil center [x y 0]
nlist = zeros(ncoil,3); % normal vector (inward) from coil center

% circular coil configuration, like head coils
alist = deg2rad(orbit)/ncoil * [0:(ncoil-1)]; % list of coil angles in radians
for ii=1:ncoil
	phi = alist(ii);
	Rad = max(nx/2 * dx, ny/2 * dy) * coil_distance;
	plist(ii,:) = Rad * [cos(phi) sin(phi) 0];
	nlist(ii,:) = -[cos(phi) sin(phi) 0];
	olist(ii,:) = [-sin(phi) cos(phi) 0]; % unit vector orthogonal to nlist
end

% object coordinates for slice z=0
x = ([1:nx] - (nx+1)/2) * dx;
y = ([1:ny] - (ny+1)/2) * dy;
z = 0;
[xx yy zz] = ndgrid(x,y,z);

smap = zeros(nx, ny, ncoil, 'single');
for ii=1:ncoil
	% rotate coordinates to correspond to coil orientation
	zr =	(xx - plist(ii,1)) .* nlist(ii,1) + ...
		(yy - plist(ii,2)) .* nlist(ii,2) + ...
		(zz - plist(ii,3)) .* nlist(ii,3);
	xr =	xx .* nlist(ii,2) - yy .* nlist(ii,1);

	if 0 % see coordinates
		im plc 1 2
		im(1, x, y, xr), xlabel x, ylabel y
		im(2, x, y, zr)
		keyboard 
	end

	[sx sy sz] = mri_smap1(xr, 0, zr, rlist(ii)); % in coil coordinates

	if 0 % see field components
		im plc 2 2
		im(1, x, y, sx), cbar
		im(2, x, y, sy), cbar
		im(3, x, y, sz), cbar
		im subplot 2
		tmp = sqrt(sx.^2 + sz.^2);
		quiver(x, y, (sx./tmp)', (sz./tmp)', 0), axis square
	end

	if flag_old
		smap(:,:,ii) = sz; % old way (wrong!) based only on smap_z

	else
		if nlist(ii,3) || olist(ii,3)
			fail 'unsupported'
		end
		% assume z component of plist and nlist are 0
		bx = sz * nlist(ii,1) + sx * olist(ii,1);
		by = sz * nlist(ii,2) + sx * olist(ii,2);
		smap(:,:,ii) = bx + 1i * by;

		if 0 % see final field components vs phase
			im subplot 4
			bb = sqrt(bx.^2 + by.^2);
			quiver(x, y, (bx./bb)', (by./bb)', 0), axis square
			im(2, x, y, angle(smap(:,:,ii))), cbar
			keyboard
		end
	end
end
smap = smap * rlist(1) / (2*pi); % trick: scale so maximum is near unity

if chat && im % show smap and array geometry in z=0 plane
	mri_sensemap_sim_show(smap, x, y, dx, dy, nlist, plist, rlist)
end


% mri_sensemap_sim_show()
function mri_sensemap_sim_show(smap, x, y, dx, dy, nlist, plist, rlist)
[nx ny ncoil] = size(smap);
nshow = min(max(ncoil,2),4);
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
	plot(0,0,'.', plist(:,1), plist(:,2), 'o')
	xdir = nlist(ii,2);
	ydir = nlist(ii,1);
	r = rlist(ii);
	plot(plist(ii,1)+r*xdir*[-1 1], plist(ii,2)+r*ydir*[1 -1], '-')
	hold off

%	ph = ir_unwrap(angle(tmp)); % trick: unwrap phase for pretty display
	ph = angle(tmp); % show raw phase (understandable with hsv colormap)
	im(ii+nshow, 'hsv', x, y, ph, [-pi pi], 'Phase'), cbar
	axis(xmax*[-1 1 -1 1]*1.1)
	xtick([-1 0 1] * nx/2 * dx)
	ytick([-1 0 1] * ny/2 * dy)
end

ssos = sqrt(sum(abs(smap), 3)); % bug! missing .^2 !
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


% mri_smap_r(r, z)
% function for testing near 0
function out = mri_smap_r(r, z)
M = 4 * r ./ ((1 + r).^2 + z.^2); % = k^2, see ellipke
[K E] = ellipke(M);
out = 2 * z ./ r .* ((1 + r).^2 + z.^2).^(-0.5) .* ...
	((1 + r.^2 + z.^2) ./ ((1 - r).^2 + z.^2) .* E - K);


% mri_smap1()
% based on grivich:00:tmf
% for a circular coil in "x-y plane" of radius a
% note that coil x-y plane is not same as object x-y plane!
function [smap_x smap_y smap_z] = mri_smap1(x, y, z, a)
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
	t0 = mri_smap_r(r0, z0);
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

%phi = atan2(y, x);
%smap_x = smap_r .* cos(phi);
%smap_y = smap_r .* sin(phi);


% mri_sensemap_sim_test1
% test mri_smap1 routine, cf Fig. 4 of grivich:00:tmf
function mri_sensemap_sim_test1
a = 1;
x = linspace(-2,2,99);
y = linspace(-2,2,97);
zlist = [0.1 0.2 0.5 1.0];
[xx yy zz] = ndgrid(x, y, zlist);
[smap_x smap_y smap_z] = mri_smap1(xx, yy, zz, a);
im('plc', 4, numel(zlist))
mri_sensemap_sim_test1_show(smap_x, x, y, 0, zlist, 'x')
mri_sensemap_sim_test1_show(smap_y, x, y, 4, zlist, 'y')
mri_sensemap_sim_test1_show(smap_z, x, y, 8, zlist, 'z')
smap_b = sqrt(smap_x.^2 + smap_y.^2);
mri_sensemap_sim_test1_show(smap_b, x, y, 12, zlist, 'b')


% mri_sensemap_sim_test1_show()
function mri_sensemap_sim_test1_show(map, x, y, offset, zlist, leg)
clim = [-20 20];
for iz = 1:numel(zlist)
	p = (iz-1) * 4;
	im(offset+iz, x, y, map(:,:,iz), leg, clim), cbar
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
end
drawnow


% mri_sensemap_sim_test
function mri_sensemap_sim_test

[smap x y] = mri_sensemap_sim('chat', 1, 'nx', 32, ...
	'rcoil', [], 'ncoil', 4, 'coil_distance', 1.2);

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

if 0 && im % see ellipke
	m = linspace(0,1,101);
	[k e] = ellipke(m);
	clf, plot(m, k, '-', m, e, '--'), legend('k', 'e')
	yaxis_pi('0 p/2 p 3*p/2')
end
