 function [phantom, params] = ellipsoid_im(ig, params, varargin)
%function [phantom, params] = ellipsoid_im(ig, params, varargin)
%|
%| Generate ellipsoids phantom image from parameters:
%| in
%|	ig		image_geom()
%|	params [ne 9]	ellipsoid parameters.  if empty, use 3d shepp-logan
%|			[x_center y_center z_center  x_radius y_radius z_radius
%|				xy_angle_degrees z_angle_degrees  amplitude]
%| option
%|	'oversample'	int	over-sampling factor (default: 1)
%|	'checkfov'	0|1 	warn if any ellipsoid is out of fov
%|	hu_scale	double	use 1000 to scale shepp-logan to HU (default: 1)
%|	'type'		char	'' (default); 'fast' new 'fast' (?) version
%|			'lowmem' uses less memory than 'fast' but a bit slower
%|
%| out
%|	phantom		[nx ny nz] image
%|
%| op ellipsoid in aspire with nint=3 is oversample=4 = 2^(3-1) here
%|
%| Copyright 2004-8-13, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
%| and Jeff Fessler, University of Michigan

if nargin == 1 && streq(ig, 'test'), ellipsoid_im_test, return, end
if nargin == 1 && streq(ig, 'profile'), ellipsoid_im_profile, return, end
if nargin == 1 && streq(ig, 'big'), ellipsoid_im_test_big, return, end
if nargin < 2, ir_usage, end

if isempty(ig)
	ig = image_geom('nx', 64, 'ny', 60, 'nz', 16, 'dx', 1, 'dz', 1);
end

if isnumeric(ig)
	[phantom, params] = ellipsoid_im_old(ig, params, varargin{:});
return
end

arg.oversample = 1;
arg.checkfov = false;
arg.type = '';
arg.showmem = false;
arg.hu_scale = 1;
arg = vararg_pair(arg, varargin);

if ~isvar('params') || isempty(params)
	params = 'shepp-logan';
end

if ischar(params)
	params = shepp_logan_3d_parameters(ig.fov/2, ig.fov/2, ig.zfov/2, params);
end
params(:,9) = params(:,9) * arg.hu_scale;

if arg.checkfov
	if ~ellipsoid_im_check_fov(ig.nx, ig.ny, ig.nz, params, ...
		ig.dx, ig.dy, ig.dz, ig.offset_x, ig.offset_y, ig.offset_z)
		warn 'ellipsoid exceeds FOV'
	end
end

args = {};
switch arg.type
case 'lowmem-slow'
	fun = @ellipsoid_im_lowmem;
	args = {'fun', @ellipsoid_im_slow};
case {'lowmem', 'lowmem-fast'}
	fun = @ellipsoid_im_lowmem;
case 'fast'
	fun = @ellipsoid_im_fast;
case {'', 'slow'}
	fun = @ellipsoid_im_slow;
otherwise
	fail('bad type %s', arg.type)
end

phantom = fun(ig.nx, ig.ny, ig.nz, params, ...
	ig.dx, ig.dy, ig.dz, ig.offset_x, ig.offset_y, ig.offset_z, ...
	arg.oversample, arg.showmem, args{:});
end % ellipsoid_im


%
% ellipsoid_im_old()
%
function [phantom, params] = ellipsoid_im_old(nx, ny, nz, params, ...
	dx, dy, dz, varargin)

if nargout > 1, warning 'units of output params not finished', end

arg.oversample = 1;
arg = vararg_pair(arg, varargin);

if ~isvar('ny') || isempty(ny), ny = nx; end
if ~isvar('nz') || isempty(nz), nz = nx; end

if ~isvar('dx') || isempty(dx), dx = 1; end
if ~isvar('dy') || isempty(dy), dy = dx; end
if ~isvar('dz') || isempty(dz), dz = dx; end

%if ~isvar('params') || isempty(params), params = 'shepp-logan'; end
if ~isvar('params') || isempty(params), fail 'not done'; end

phantom = ellipsoid_im_slow(nx, ny, nz, params, ...
	dx, dy, dz, 0, 0, 0, arg.oversample, false, false);
end % ellipsoid_im_old()


%
% ellipsoid_im_slow()
%
function phantom = ellipsoid_im_slow(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, over, showmem, varargin)
if length(varargin), error 'bug', end

if size(params,2) ~= 9
	fail('bad ellipse parameter vector size')
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
	cz = par(3);	rz = par(6);

	azim = deg2rad(par(7));
	polar = deg2rad(par(8));
	[xr yr zr] = rot3(xx-cx, yy-cy, zz-cz, azim, polar);

	tmp = (xr / rx).^2 + (yr / ry).^2 + (zr / rz).^2 <= 1;
	phantom = phantom + par(9) * tmp;
end

if showmem, jf whos, end
phantom = downsample3(phantom, over);
end % ellipsoid_im_slow()


%
% ellipsoid_im_fast()
%
function phantom = ellipsoid_im_fast(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, ...
	over, showmem, varargin)
if length(varargin), fail('bug'), end

if size(params,2) ~= 9
	fail('bad ellipse parameter vector size')
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
[xx, yy, zz] = ndgrid(xx, yy, zz);

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
for ip = 1:np;
	ticker(mfilename, ip, np)

	par = params(ip, :);
	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	rz = par(6);

	azim = deg2rad(par(7));
	polar = deg2rad(par(8));

	xs = xx - cx; % shift per center
	ys = yy - cy;
	zs = zz - cz;

	[xr yr zr] = rot3(xs, ys, zs, azim, polar);
	if over == 1
		vi = (xr / rx).^2 + (yr / ry).^2 + (zr / rz).^2 <= 1;
		phantom = phantom + par(9) * single(vi);
	continue
	end

	% check all 8 corners of the cube that bounds all possible
	% 3D rotations of the voxel
	hh = sqrt(hx^2 + hy^2 + hz^2);
	vi = true; % voxels entirely inside the ellipsoid
%	vo = true; % voxels entirely outside the ellipsoid
	for iz = [-1 1]
	for iy = [-1 1]
	for ix = [-1 1]
%		pr [ix iy iz]
		xo = xr + ix * hh;
		yo = yr + iy * hh;
		zo = zr + iz * hh;
		vi = vi & ((xo / rx).^2 + (yo / ry).^2 + (zo / rz).^2 < 1);
%		vo = vo & ((xo / rx).^2 + (yo / ry).^2 + (zo / rz).^2 > 1);
	end
	end
	end

%	vo_tmp = vo; 
	vo = false; % todo: for now, "outside" test is failing
	if ip == 1, warn 'todo: must debug this', end

	if any(vi(:) & vo(:)), fail 'bug', end

if 0
	% coordinates of "outer" corner of each voxel, relative to ellipsoid center
	sign_p = @(x) (x >= 0) * sqrt(3); % modified sign()
	xo = xr + sign_p(xr) * hx;
	yo = yr + sign_p(yr) * hy;
	zo = zr + sign_p(zr) * hz;

	% voxels that are entirely inside the ellipsoid
	vi = (xo / rx).^2 + (yo / ry).^2 + (zo / rz).^2 <= 1;
end

if 0
	% coordinates of "inner" corner of each pixel, relative to ellipse center
	sign_n = @(x) (x > 0) * sqrt(3); % modified sign()
	xi = xr - sign_n(xs) * hx;
	yi = yr - sign_n(ys) * hy;
	zi = zr - sign_n(zs) * hz;

	% voxels that are entirely outside the ellipsoid
	vo = (max(abs(xi),0) / rx).^2 + (max(abs(yi),0) / ry).^2 ...
		+ (max(abs(zi),0) / rz).^2 >= 1;
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

	[xr, yr, zr] = rot3(x, y, z, azim, polar);
	in = (xr / rx).^2 + (yr / ry).^2 + (zr / rz).^2 <= 1;
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

	phantom = phantom + par(9) * gray;
end

if showmem, jf whos, end

end % ellipsoid_im_fast()


%
% ellipsoid_im_lowmem()
% this version does "one slice at a time" to reduce memory
%
function phantom = ellipsoid_im_lowmem(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, over, showmem, varargin)
arg.fun = @ellipsoid_im_fast;
arg = vararg_pair(arg, varargin);

phantom = zeros(nx, ny, nz, 'single');
for iz=1:nz
	offset_z_new = (nz-1)/2 + offset_z - (iz-1);
	phantom(:,:,iz) = arg.fun(nx, ny, 1, params, ...
		dx, dy, dz, offset_x, offset_y, offset_z_new, over, ...
		showmem && iz == 1);
end

end % ellipsoid_im_lowmem()



%
% ellipsoid_im_check_fov()
%
function ok = ellipsoid_im_check_fov(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z)

if size(params,2) ~= 9
	fail('bad ellipse parameter vector size')
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
	cz = par(3);	rz = par(6);

	if cx + rx > xmax || cx - rx < xmin
		warn('fov: x range %g %g, cx=%g rx=%g', xmin, xmax, cx, rx)
		ok = false;
	end
	if cy + ry > ymax || cy - ry < ymin
		warn('fov: y range %g %g, cy=%g ry=%g', ymin, ymax, cy, ry)
		ok = false;
	end
	if cz + rz > zmax || cz - rz < zmin
		warn('fov: z range %g %g, cz=%g rz=%g', zmin, zmax, cz, rz)
		ok = false;
	end
end

end % ellipsoid_im_check_fov


%
% rot3()
%
function [xr, yr, zr] = rot3(x, y, z, azim, polar)
if polar, error 'z (polar) rotation not done', end
xr =  cos(azim) * x + sin(azim) * y;
yr = -sin(azim) * x + cos(azim) * y;
zr = z;
end % rot3()


%
% shepp_logan_3d_parameters()
% most of these values are unitless "fractions of field of view"
%
function params = shepp_logan_3d_parameters(xfov, yfov, zfov, ptype)

% parameters from Kak and Slaney text, p. 102, which seem to have typos!
ekak = [...
	0	0	0	0.69	0.92	0.9	0	2.0;
	0	0	0	0.6624	0.874	0.88	0	-0.98;
	-0.22	0	-0.25	0.41	0.16	0.21	108	-0.02;
	0.22	0	-0.25	0.31	0.11	0.22	72	-0.02;
	0	0.1	-0.25	0.046	0.046	0.046	0	0.02; % same?
	0	0.1	-0.25	0.046	0.046	0.046	0	0.02; % same?
	-0.8	-0.65	-0.25	0.046	0.023	0.02	0	0.01;
	0.06	-0.065	-0.25	0.046	0.023	0.02	90	0.01;
	0.06	-0.105	0.625	0.56	0.04	0.1	90	0.02;
	0	0.1	-0.625	0.056	0.056	0.1	0	-0.02];

% the following parameters came from leizhu@stanford.edu
% who says that the Kak&Slaney values are incorrect
% fix: i haven't had time to look into this in detail
% yu:05:ads cites shepp:74:tfr 

%	x	y	z	rx	ry	rz	angle	density
ezhu = [...
	0	0	0	0.69	0.92	0.9	0	2.0;
	0	-0.0184	0	0.6624	0.874	0.88	0	-0.98;
	-0.22	0	-0.25	0.41	0.16	0.21	-72	-0.02;
	0.22	0	-0.25	0.31	0.11	0.22	72	-0.02;
	0	0.35	-0.25	0.21	0.25	0.35	0	0.01;
	0	0.1	-0.25	0.046	0.046	0.046	0	0.01;
	-0.08	-0.605	-0.25	0.046	0.023	0.02	0	0.01;
	0	-0.1	-0.25	0.046	0.046	0.046	0	0.01;
	0	-0.605	-0.25	0.023	0.023	0.023	0	0.01;
	0.06	-0.605	-0.25	0.046	0.023	0.02	-90	0.01;
	0.06	-0.105	0.0625	0.056	0.04	0.1	-90	0.02;
	0	0.1	0.625	0.056	0.056	0.1	0	-0.02];

% and here are parameters from the "phantom3d.m" in matlab central
% http://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom 
% by Matthias Schabel matlab@schabel-family.org
% citing p199-200 of peter toft thesis: http://petertoft.dk/PhD/
% but that thesis has only 2d phantom!
%
% e(:,1) = [1 -.98 -.02 -.02 .01 .01 .01 .01 .01 .01];
%
% 1:	A	additive intensity value of the ellipsoid
% 2:	a	length of the x semi-axis of the ellipsoid 
% 3:	b	length of the y semi-axis of the ellipsoid
% 4:	c	length of the z semi-axis of the ellipsoid
% 5:	x0	x-coordinate of the center of the ellipsoid
% 6:	y0	y-coordinate of the center of the ellipsoid
% 7:	z0	z-coordinate of the center of the ellipsoid
% 8:	phi	Euler angle (in degrees) (rotation about z-axis)
% 9:	theta	Euler angle (in degrees) (rotation about x-axis)
% 10:	psi	Euler angle (in degrees) (rotation about z-axis)
%
% For purposes of generating the phantom, the domains for the x-, y-, and 
% z-axes span [-1,1].  Columns 2 through 7 must be specified in terms of
% this range.
%
%         A     a    b    c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e3d =  [  1 .6900 .920 .810      0       0       0      0      0      0
        -.8 .6624 .874 .780      0  -.0184       0      0      0      0
        -.2 .1100 .310 .220    .22       0       0    -18      0     10
        -.2 .1600 .410 .280   -.22       0       0     18      0     10
         .1 .2100 .250 .410      0     .35    -.15      0      0      0
         .1 .0460 .046 .050      0      .1     .25      0      0      0
         .1 .0460 .046 .050      0     -.1     .25      0      0      0
         .1 .0460 .023 .050   -.08   -.605       0      0      0      0
         .1 .0230 .023 .020      0   -.606       0      0      0      0
         .1 .0230 .046 .020    .06   -.605       0      0      0      0 ];

switch ptype
case {'shepp-logan', 'shepp-logan-zhu', 'zhu', ''}
	params = ezhu;
case {'shepp-logan-kak', 'kak'}
	params = ekak;
case {'shepp-logan-e3d', 'e3d'}
	params = e3d(:, [5:7 2:4 8 1]); % x y z rx ry rz angle density
otherwise
	fail('unknown parameter type %s', ptype)
end

params(:,[1 4]) = params(:,[1 4]) * xfov;
params(:,[2 5]) = params(:,[2 5]) * yfov;
params(:,[3 6]) = params(:,[3 6]) * zfov;
params(:,9) = params(:,8);
params(:,8) = 0; % z rotation
end % shepp_logan_3d_parameters()


%
% ellipsoid_im_profile()
%
function ellipsoid_im_profile
ig = image_geom('nx', 2^6, 'ny', 2^6-2, 'nz', 2^5, 'fov', 240, 'dz', 1);
profile on
phantom = ellipsoid_im(ig, [], 'oversample', 2, 'type', '');
phantom = ellipsoid_im(ig, [], 'oversample', 2, 'type', 'fast');
phantom = ellipsoid_im(ig, [], 'oversample', 2, 'type', 'lowmem');
profile off
profile report

end % ellipsoid_im_profile()


%
% ellipsoid_im_test_big()
% compare timings for big size
%
function ellipsoid_im_test_big

ellipsoid_im([], []); % warm up
ig = image_geom('nx', 2^8, 'ny', 2^8-2^4, 'nz', 2^7, ...
	'fov', 400, 'dz', 1.1, 'down', 1);
types = {'fast', 'lowmem'};
for ii=1:length(types)
	type = types{ii};
	cpu etic
	mat{ii} = ellipsoid_im(ig, [], 'oversample', 3, ...
		'type', type, 'showmem', true);
	cpu('etoc', sprintf('type=%s', type))
	if ii > 1
		max_percent_diff(mat{1}, mat{ii})
	end
end
end % ellipsoid_im_test_big


%
% ellipsoid_im_test()
%
function ellipsoid_im_test

if 1
	ig = image_geom('nx', 512, 'fov', 500, 'nz', 64, 'dz', 0.625, 'down', 8);
	x = ellipsoid_im(ig, [], 'hu_scale', 1000);
	im clf, im(x, [900 1100]), cbar
prompt
end

ig = image_geom('nx', 2^5, 'ny', 2^5-2, 'nz', 15, 'fov', 240, ...
	'dz', -6); % negative dz to match aspire
im pl 2 2
if 1
	over = 2^0;
%	par = [30 20 10, 50 40 30, 20 0 100];
	par = [];
	args = {ig, par, 'oversample', over};
	[phantom par] = ellipsoid_im(args{:}, 'hu_scale', -5000);

	over = 2^2;
	par = par(3,:);
	par(7) = 0;
	args = {ig, par, 'oversample', over};
	[phantom par] = ellipsoid_im(args{:});
%	im(1, phantom, 'Shepp Logan', [0.9 1.1]), cbar
end

if 1 % compare different methods
	types = {'fast', 'lowmem', 'lowmem-slow'};
	for ii=1:length(types)
		tmp = ellipsoid_im(args{:}, 'type', types{ii});
		max_percent_diff(phantom, tmp)
		im(1, phantom)
		im(2, tmp)
		im(3, phantom-tmp)
%		pr unique(phantom)
%		pr unique(tmp)
	end
%return % todo: fix bug
	im(1, phantom, 'Shepp Logan', [0.9 1.1]), cbar
end

% compare to aspire
if ~has_aspire, return, end

par = [30 20 10, 50 40 30, 20 0 100];

dir = test_dir;
file = [dir '/t.fld'];
com = 'echo y | op ellipsoid %s %d %d %d  %g %g %g  %g %g %g %g %g %d %d';
pix = [ig.dx -ig.dy -ig.dz ig.dx ig.dy -ig.dz 1 1 1];
com = sprintf(com, file, ig.nx, ig.ny, ig.nz, par ./ pix, log2(over)+1);
os_run(com)
asp = fld_read(file);
im(4, ig.x, ig.y, asp, 'aspire'), cbar

mat = ellipsoid_im(ig, par, 'oversample', over, 'type', 'fast'); % todo!

if 0 % volume
	pr 4/3 * pi * prod(par(4:6)) * par(9)
	pr sum(mat(:)) * abs(ig.dx * ig.dy * ig.dz);
end

t = sprintf('mat, z: %g to %g', ig.z([1 end]));
im(2, ig.x, ig.y, mat, t), cbar
im(3, ig.x, ig.y, mat-asp, 'mat-aspire'), cbar
max_percent_diff(mat, asp) % 25% different it seems
%pr unique(mat)
%pr unique(asp)

if 1 % check centroid
	[xx yy zz] = ndgrid(ig.x, ig.y, ig.z);
	t = [sum(xx(:) .* mat(:)) sum(yy(:) .* mat(:)) ...
		sum(zz(:) .* mat(:))] / sum(mat(:));
	if any(abs(t - par(1:3)) > 0.02)
		pr par(1:3), pr t
		warn 'bad centroid'
	end
end
end % ellipsoid_im_test()
