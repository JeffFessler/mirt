  function [phantom, params] = cylinder_im(ig, params, varargin)
%|function [phantom, params] = cylinder_im(ig, params, varargin)
%|
%| Generate (elliptical) cylinder phantom image from parameters.
%| in
%|	ig		image_geom()
%|	params [ne 8]	cylinder parameters.  if empty, use 3d Defrise
%|			[x_center y_center z_center  x_radius y_radius z_length
%|				xy_angle_degrees  amplitude]
%| option
%|	'oversample'	int	over-sampling factor (default: 1)
%|	'checkfov'	0|1 	warn if any cylinder is out of fov
%|	'type'		char	'' (default); 'fast' new 'fast' (?) version
%|			'lowmem' uses less memory than 'fast' but a bit slower
%|
%| out
%|	phantom		[nx ny nz] image
%|
%| todo: could make it exact in z direction!
%|
%| Copyright 2013-06-10, Jeff Fessler, University of Michigan

if nargin == 1 && streq(ig, 'test'), cylinder_im_test, return, end
if nargin == 1 && streq(ig, 'profile'), cylinder_im_profile, return, end
if nargin == 1 && streq(ig, 'big'), cylinder_im_test_big, return, end
if nargin < 2, help(mfilename), error(mfilename), end

if isempty(ig)
	ig = image_geom('nx', 64, 'ny', 60, 'nz', 16, 'dx', 1, 'dz', 1);
end

arg.oversample = 1;
arg.checkfov = false;
arg.type = '';
arg.showmem = false;
arg = vararg_pair(arg, varargin);

if ~isvar('params') || isempty(params)
	params = 'defrise';
end

if ischar(params)
	params = cylinder_im_parameters(ig.fov/2, ig.fov/2, ig.zfov/2, params);
end

if arg.checkfov
	if ~cylinder_im_check_fov(ig.nx, ig.ny, ig.nz, params, ...
		ig.dx, ig.dy, ig.dz, ig.offset_x, ig.offset_y, ig.offset_z)
		warn 'cylinder exceeds FOV'
	end
end

args = {};
switch arg.type
case 'lowmem-slow'
	fun = @cylinder_im_lowmem;
	args = {'fun', @cylinder_im_slow};
case {'lowmem', 'lowmem-fast'}
	fun = @cylinder_im_lowmem;
case 'fast'
	fun = @cylinder_im_fast;
case {'', 'slow'}
	fun = @cylinder_im_slow;
otherwise
	fail('bad type %s', arg.type)
end

phantom = fun(ig.nx, ig.ny, ig.nz, params, ...
	ig.dx, ig.dy, ig.dz, ig.offset_x, ig.offset_y, ig.offset_z, ...
	arg.oversample, arg.showmem, args{:});
end % cylinder_im




%
% cylinder_im_slow()
%
function phantom = cylinder_im_slow(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, over, showmem, varargin)
if length(varargin), error 'bug', end

if size(params,2) ~= 8
	error 'bad cylinder parameter vector size'
end

phantom = zeros(nx*over, ny*over, nz*over, 'single');

wx = (nx*over-1)/2 + offset_x*over;
wy = (ny*over-1)/2 + offset_y*over;
wz = (nz*over-1)/2 + offset_z*over;
xx = ((0:nx*over-1) - wx) * dx / over;
yy = ((0:ny*over-1) - wy) * dy / over;
zz = ((0:nz*over-1) - wz) * dz / over;
xmax = max(xx); xmin = min(xx);
ymax = max(yy); ymin = min(yy);
zmax = max(zz); zmin = min(zz);
[xx yy zz] = ndgrid(xx, yy, zz);

ticker reset
np = nrow(params);
for ip = 1:np;
	ticker(mfilename, ip, np)

	par = params(ip, :);
	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	zh = par(6) / 2;

	azim = deg2rad(par(7));
	polar = 0;
	[xr yr zr] = rot3(xx-cx, yy-cy, zz-cz, azim, polar);

	tmp = ((xr / rx).^2 + (yr / ry).^2 <= 1) & abs(zr / zh) <= 1;
	phantom = phantom + par(8) * tmp;
end

if showmem, jf whos, end
phantom = downsample3(phantom, over);
end % cylinder_im_slow()


%
% cylinder_im_fast()
%
function phantom = cylinder_im_fast(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, ...
	over, showmem, varargin)
if length(varargin), error 'bug', end

if size(params,2) ~= 8
	error 'bad cylinder parameter vector size'
end

phantom = zeros(nx, ny, nz, 'single');

wx = (nx-1)/2 + offset_x;
wy = (ny-1)/2 + offset_y;
wz = (nz-1)/2 + offset_z;
xx = ((0:nx-1) - wx) * dx;
yy = ((0:ny-1) - wy) * dy;
zz = ((0:nz-1) - wz) * dz;
xmax = max(xx); xmin = min(xx);
ymax = max(yy); ymin = min(yy);
zmax = max(zz); zmin = min(zz);
[xx yy zz] = ndgrid(xx, yy, zz);

if over > 1
	tmp = ((1:over) - (over+1)/2) / over;
	[xf yf zf] = ndgrid(tmp*dx, tmp*dy, tmp*dz);
	xf = xf(:)';
	yf = yf(:)';
	zf = zf(:)';

	hx = abs(dx) / 2;
	hy = abs(dy) / 2;
	hz = abs(dz) / 2;
end

ticker reset
np = nrow(params);
for ip = 1:np
	ticker(mfilename, ip, np)

	par = params(ip, :);
	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	zh = par(6) / 2;

	azim = deg2rad(par(7));
	polar = 0;

	xs = xx - cx; % shift per center
	ys = yy - cy;
	zs = zz - cz;

	[xr yr zr] = rot3(xs, ys, zs, azim, polar);
	if over == 1
		vi = ((xr / rx).^2 + (yr / ry).^2 <= 1) & abs(zr / zh) <= 1;
		phantom = phantom + par(8) * single(vi);
	continue
	end

	% check all 8 corners of the cube that bounds all possible
	% 3D rotations of the voxel
	hh = sqrt(hx^2 + hy^2 + hz^2);
	vi = true; % voxels entirely inside the cylinder
	for iz = [-1 1]
	for iy = [-1 1]
	for ix = [-1 1]
%		pr [ix iy iz]
		xo = xr + ix * hh;
		yo = yr + iy * hh;
		zo = zr + iz * hh;
		vi = vi & ((xo / rx).^2 + (yo / ry).^2 <= 1) & abs(zo / zh) <= 1;
	end
	end
	end

	vo = false; % todo: for now, "outside" test is failing
	if ip == 1, warn 'todo: must debug this', end

	if any(vi(:) & vo(:)), fail 'bug', end

if 0
	% coordinates of "outer" corner of each voxel, relative to cylinder center
	sign_p = @(x) (x >= 0) * sqrt(3); % modified sign()
	xo = xr + sign_p(xr) * hx;
	yo = yr + sign_p(yr) * hy;
	zo = zr + sign_p(zr) * hz;

	% voxels that are entirely inside the cylinder
	vi = ((xo / rx).^2 + (yo / ry).^2 <= 1) & abs(zo / zh) <= 1;
end

	% subsampling for edge voxels
	edge = ~vi & ~vo;
%	edge = true(size(edge)); % for testing
	if showmem
		printm('edge fraction %g = %d / %d', ...
			sum(edge(:)) / prod(size(edge)), ...
			sum(edge(:)),prod(size(edge)))
	end
	x = xx(edge) - cx;
	y = yy(edge) - cy;
	z = zz(edge) - cz;
	x = outer_sum(x, xf);
	y = outer_sum(y, yf);
	z = outer_sum(z, zf);

	[xr yr zr] = rot3(x, y, z, azim, polar);
	in = ((xr / rx).^2 + (yr / ry).^2 <= 1) & abs(zr / zh) <= 1;
	tmp = mean(in, 2);

	gray = single(vi);
	gray(edge) = tmp;

	if 0 % todo: help debug
		ee = (gray > 0) & (gray < 1);
		im(ee)
		im(vi)
	prompt
		%minmax(tmp)
		%clf, im pl 1 2
		%im(tmp)
	end

	phantom = phantom + par(8) * gray;
end

if showmem, jf whos, end

end % cylinder_im_fast()


%
% cylinder_im_lowmem()
% this version does "one slice at a time" to reduce memory
%
function phantom = cylinder_im_lowmem(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, over, showmem, varargin)
arg.fun = @cylinder_im_fast;
arg = vararg_pair(arg, varargin);

phantom = zeros(nx, ny, nz, 'single');
for iz=1:nz
	offset_z_new = (nz-1)/2 + offset_z - (iz-1);
	phantom(:,:,iz) = arg.fun(nx, ny, 1, params, ...
		dx, dy, dz, offset_x, offset_y, offset_z_new, over, ...
		showmem && iz == 1);
end

end % cylinder_im_lowmem()



%
% cylinder_im_check_fov()
%
function ok = cylinder_im_check_fov(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z)

if size(params,2) ~= 8
	error 'bad cylinder parameter vector size'
end

wx = (nx-1)/2 + offset_x;
wy = (ny-1)/2 + offset_y;
wz = (nz-1)/2 + offset_z;
xx = ((0:nx-1) - wx) * dx;
yy = ((0:ny-1) - wy) * dy;
zz = ((0:nz-1) - wz) * dz;
xmax = max(xx); xmin = min(xx);
ymax = max(yy); ymin = min(yy);
zmax = max(zz); zmin = min(zz);

ok = true;
for ip = 1:nrow(params)
	par = params(ip, :);
	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	zh = par(6) / 2;

	if cx + rx > xmax || cx - rx < xmin
		warn('fov: x range %g %g, cx=%g rx=%g', xmin, xmax, cx, rx)
		ok = false;
	end
	if cy + ry > ymax || cy - ry < ymin
		warn('fov: y range %g %g, cy=%g ry=%g', ymin, ymax, cy, ry)
		ok = false;
	end
	if cz + zh > zmax || cz - zh < zmin
		warn('fov: z range %g %g, cz=%g zh=%g', zmin, zmax, cz, zh)
		ok = false;
	end
end

end % cylinder_im_check_fov


%
% rot3()
%
function [xr, yr, zr] = rot3(x, y, z, azim, polar)
if polar, error 'z (polar) rotation not done', end
xr =  cos(azim) * x + sin(azim) * y;
yr = -sin(azim) * x + cos(azim) * y;
zr = z;
end % rot3()



% cylinder_im_parameters()
function params = cylinder_im_parameters(xfov, yfov, zfov, ptype)
switch ptype
case 'defrise'
	params = [ ...
		[0.2 0.1 -0.75	0.5 0.6 0.2	20	1000];
		[0.2 0.1 -0.25	0.5 0.6 0.2	20	1000];
		[0.2 0.1 0.25	0.5 0.6 0.2	20	1000];
		[0.2 0.1 0.75	0.5 0.6 0.2	20	1000];
	];

otherwise
	fail 'unknown'
end

params(:,[1 4]) = params(:,[1 4]) * xfov;
params(:,[2 5]) = params(:,[2 5]) * yfov;
params(:,[3 6]) = params(:,[3 6]) * zfov;
end % cylinder_im_parameters()


%
% cylinder_im_profile()
%
function cylinder_im_profile
ig = image_geom('nx', 2^6, 'ny', 2^6-2, 'nz', 2^5, 'fov', 240, 'dz', 1);
profile on
phantom = cylinder_im(ig, [], 'oversample', 2, 'type', '');
phantom = cylinder_im(ig, [], 'oversample', 2, 'type', 'fast');
phantom = cylinder_im(ig, [], 'oversample', 2, 'type', 'lowmem');
profile off
profile report

end % cylinder_im_profile()


%
% cylinder_im_test_big()
% compare timings for big size
%
function cylinder_im_test_big

cylinder_im([], []); % warm up
ig = image_geom('nx', 2^8, 'ny', 2^8-2^4, 'nz', 2^7, ...
	'fov', 400, 'dz', 1.1, 'down', 1);
types = {'fast', 'lowmem'};
for ii=1:length(types)
	type = types{ii};
	cpu etic
	mat{ii} = cylinder_im(ig, [], 'oversample', 3, ...
		'type', type, 'showmem', true);
	cpu('etoc', sprintf('type=%s', type))
	if ii > 1
		max_percent_diff(mat{1}, mat{ii})
	end
end
end % cylinder_im_test_big


%
% cylinder_im_test()
%
function cylinder_im_test

if 1
	ig = image_geom('nx', 512, 'fov', 500, 'nz', 64*8, 'dz', 0.625, 'down', 8);
	x = cylinder_im(ig, [], 'oversample', 1);
	im clf, im(x), cbar
%	im clf, im(permute(x, [1 3 2])), cbar
prompt
end

if 1
	ig = image_geom('nx', 2^5, 'ny', 2^5-2, 'nz', 15, 'fov', 240, ...
		'dz', -6); % negative dz to match aspire
	args = {ig, 'defrise', 'oversample', 2^1};
	phantom = cylinder_im(args{:});

	im plc 2 4
	im(1, phantom, 'Defrise')
end

if 1 % compare different methods
	types = {'fast', 'lowmem', 'lowmem-slow'};
	for ii=1:numel(types)
		tmp = cylinder_im(args{:}, 'type', types{ii});
		max_percent_diff(phantom, tmp)
		im(ii+1, tmp)
		im(ii+5, phantom-tmp)
	end
end

end % cylinder_im_test()
