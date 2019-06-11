 function [phantom, params] = cuboid_im(ig, params, varargin)
%function [phantom, params] = cuboid_im(ig, params, varargin)
%|
%| generate cuboid phantom image from parameters:
%|	[x_center y_center z_center  x_diameter y_diameter z_diameter
%|		xy_angle_degrees z_angle_degrees amplitude]
%| in
%|	ig			image_geom()
%|	params		[9 N]	cuboid parameters.  if empty, use default
%|				Note!! "diameter" not "radius"
%| option
%|	'oversample'	[]	over-sampling factor
%|	'type'		char	'sample' use samples
%|				'lowmem1' one slice per time
%|				'exact' perfect partial volume if angle* = 0
%|				'' (default): 'exact' if angle* = 0 else 'samples'
%| out
%|	phantom		[nx ny nz] image
%|
%| Yong Long, 2008-08-28, University of Michigan adapted from ellipsoid_im()
%| 2013-05-24 Jeff Fessler
%| 2013-06-28 Jeff Fessler, added 'exact' option for non-rotated cuboids

if nargin == 1 && streq(ig, 'test'), cuboid_im_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.oversample = 1;
arg.type = '';
arg = vararg_pair(arg, varargin);

if ~isvar('params') || isempty(params)
	params = [0 0 0 abs([2*ig.dx 2*ig.dy 2*ig.dz]) 0 0 1];
end

if size(params,2) ~= 9
	fail 'bad cuboid parameter vector size'
end

if isempty(arg.type)
	if all(col(params(:,[7 8]) == 0))
		arg.type = 'exact';
	else
		arg.type = 'sample';
	end
end

if any(col(params(:,4:6) < 0))
	fail 'cuboid "diameters" must be nonnegative'
end

switch arg.type
case 'exact'
	fun = @cuboid_im_do_exact;
	if arg.oversample ~= 1 % pointless for exact version
		warn('ignoring oversample=%g for "exact"', arg.oversample)
	end
case 'lowmem1'
	fun = @cuboid_im_do_lowmem1;
case 'sample'
	fun = @cuboid_im_do;
otherwise
	fail('bad type %s', arg.type)
end

[phantom params] = fun(ig.nx, ig.ny, ig.nz, params, ...
	ig.dx, ig.dy, ig.dz, ig.offset_x, ig.offset_y, ig.offset_z, ...
	arg.oversample);

end % cuboid_im


% cuboid_im_do()
function [phantom params] = cuboid_im_do(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, ...
	over)

phantom = zeros(nx*over, ny*over, nz*over);

wx = (nx*over-1)/2 + offset_x*over;
wy = (ny*over-1)/2 + offset_y*over;
wz = (nz*over-1)/2 + offset_z*over;
xx = ([0:nx*over-1] - wx) * dx / over;
yy = ([0:ny*over-1] - wy) * dy / over;
zz = ([0:nz*over-1] - wz) * dz / over;
[xx yy zz] = ndgrid(xx, yy, zz);

ticker reset
ne = nrow(params);
for ie = 1:ne;
	ticker(mfilename, ie, ne)

	par = params(ie, :);
	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	rz = par(6);

	theta = deg2rad(par(7));
	phi = deg2rad(par(8));
	if phi, error 'z rotation not done', end
	x = cos(theta) * (xx-cx) + sin(theta) * (yy-cy);
	y = -sin(theta) * (xx-cx) + cos(theta) * (yy-cy);
	z = zz - cz;
	% rx, ry and rz are "diameters" not "radius"
	tmp = abs(x/rx) <= 1/2 & abs(y/ry) <= 1/2 & abs(z/rz) <= 1/2;

	phantom = phantom + par(9) * tmp;
end

phantom = downsample3(phantom, over);
end % cuboid_im_do()


% cuboid_im_do_exact()
function [phantom params] = cuboid_im_do_exact(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, ...
	over_this_is_ignored)

phantom = zeros(nx, ny, nz);

wx = (nx-1)/2 + offset_x;
wy = (ny-1)/2 + offset_y;
wz = (nz-1)/2 + offset_z;
xx = ([0:nx-1] - wx) * dx;
yy = ([0:ny-1] - wy) * dy;
zz = ([0:nz-1] - wz) * dz;
[xx yy zz] = ndgrid(xx, yy, zz);

% length of intersection of [a,b] with [c,d]
fun1 = @(a,b,c,d) max(min(d,b) - max(a,c), 0);
%fun2 = @(a,b,c,d) fun1(min(a,b), max(a,b), min(c,d), max(c,d));
fun2 = @(x,h) fun1(x-0.5, x+0.5, -h, h);

adx = abs(dx);
ady = abs(dy);
adz = abs(dz);

ticker reset
ne = nrow(params);
for ie = 1:ne;
	ticker(mfilename, ie, ne)

	par = params(ie, :);
	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	rz = par(6);

	if par(7) || par(8), error 'rotation not done', end

	% rx, ry and rz are "diameters" not "radius"
	x = (xx - cx) / adx; hx2 = rx / adx / 2;
	y = (yy - cy) / ady; hy2 = ry / ady / 2;
	z = (zz - cz) / adz; hz2 = rz / adz / 2;

	fx = fun2(x, hx2);
	fy = fun2(y, hy2);
	fz = fun2(z, hz2);
	tmp = (fx .* fy .* fz);
	phantom = phantom + par(9) * tmp;
end

end % cuboid_im_do_exact()


% cuboid_im_do_lowmem1()
% this version does "one slice at a time" to reduce memory
function [phantom params] = cuboid_im_do_lowmem1(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, over)

phantom = zeros(nx, ny, nz, 'single');
for iz=1:nz
	offset_z_new = (nz-1)/2 + offset_z - (iz-1);
	phantom(:,:,iz) = cuboid_im_do(nx, ny, 1, params, ...
		dx, dy, dz, offset_x, offset_y, offset_z_new, over);
end

end % cuboid_im_lowmem1()


% cuboid_im_test()
function cuboid_im_test
ig = image_geom('nx', 2^3, 'ny', 2^3, 'nz', 2^4, 'fov', 240);

if 1 % check exact for unaligned case
	diam = abs([2*ig.dx 2.7*ig.dy 3.2*ig.dz]);
	params = [1.4 -0.5 1 diam 0 0 1];
	phantom_exact = cuboid_im(ig, params, 'type', 'exact');
	im(phantom_exact, 'cuboid'), cbar
	dxyz = abs(ig.dx * ig.dy * ig.dz);
	vol_true = prod(diam); % exact volume
	vol_phantom = sum(phantom_exact(:)) * dxyz;
	equivs(vol_phantom, vol_true)

	phantom2 = cuboid_im(ig, params, 'type', 'sample', 'oversample', 2);
	vol_phantom2 = sum(phantom2(:)) * dxyz;
	max_percent_diff(vol_phantom2, vol_true)
	pr '[vol_phantom2 vol_true vol_phantom]'

	im plc 1 3
	im(1, phantom_exact), cbar
	im(2, phantom2), cbar
	im(3, phantom2 - phantom_exact), cbar
end

if 1 % aligned case
	phantom0 = cuboid_im(ig, [], 'type', 'exact');
	phantom2 = cuboid_im(ig, [], 'oversample', 2, 'type', 'sample');
	phantom4 = cuboid_im(ig, [], 'oversample', 2, 'type', 'lowmem1');

	vol0 = sum(phantom0(:));
	vol2 = sum(phantom2(:));
	vol4 = sum(phantom4(:));
	jf_equal(2^3, vol0)
	jf_equal(2^3, vol2)
	jf_equal(2^3, vol4)
%	im(phantom0, 'default cuboid'), cbar
end

end % cuboid_im_test()
