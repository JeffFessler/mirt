 function [phantom, params] = ellipse_im(ig, params, varargin)
%| Generate ellipse phantom image from parameters.
%function [phantom, params] = ellipse_im(ig, params, options)
%|
%| Parameters: [x_center y_center x_radius y_radius angle_degrees amplitude]
%|
%| in
%|	ig		strum	image_geom() object
%|	params		[ne 6]	ellipse parameters, if empty, use shepp-logan
%|				if 'shepplogan-emis' then emission version.
%|				'shepplogan-brainweb' for brainweb-like version.
%|
%| options
%|	'rot'		float	rotate ellipses by this amount [degrees]
%|	'oversample'	int	oversampling factor, for grayscale boundaries
%|	hu_scale	float	use 1000 to scale shepp-logan to HU (default: 1)
%|	'type'		char	'' | 'slow' (default); 'fast' for new method
%|
%| out
%|	phantom		[nx ny]	image
%|
%| note: op ellipse in aspire with nsub=3 is oversample=4 = 2^(3-1) here
%|
%| trick: ellipse_im(256) makes a 256^2 shepplogan-emis image.
%|
%| Copyright 2006-2-2, Jeff Fessler, University of Michigan

if nargin == 1 && streq(ig, 'test'), ellipse_im_test, return, end
if nargin == 1 && streq(ig, 'profile'), ellipse_im_profile, return, end
if nargin < 1, ir_usage, end

if isnumeric(ig)
	switch numel(ig)
	case 1
		ig = image_geom('nx', ig, 'dx', 1);
	case 2
		ig = image_geom('nx', ig(1), 'ny', ig(2), 'dx', 1);
	otherwise
		fail('numeric ig must be [nx] or [nx ny]')
	end
	if ~isvar('params') || isempty(params)
		params = 'shepplogan-emis';
	end
	[phantom, params] = ellipse_im(ig, params, varargin{:});
%	warn 'obsolete usage; please examine "help ellipse_im"'
%	[phantom, params] = ellipse_im_old(ig, params, varargin{:});
return
end

arg.rot = 0;
arg.oversample = 1;
arg.replace = 0;
arg.hu_scale = 1;
arg.fov = [];
arg.type = '';
arg = vararg_pair(arg, varargin);
if isempty(arg.fov), arg.fov = ig.fov; end

if ~isvar('params') || isempty(params)
	params = shepp_logan_parameters(arg.fov, arg.fov);
elseif streq(params, 'shepplogan-brainweb')
	params = shepp_logan_parameters(arg.fov, arg.fov);
	params(:,6) = [1 0 2 3 4 5 6 7 8 9]'; % brainweb uses index 1-10
	arg.replace = 0;
elseif streq(params, 'shepplogan-emis')
	params = shepp_logan_parameters(arg.fov, arg.fov);
%	params(:,6) = [1 2 0 4 5 6 7 8 2 2]';
%	arg.replace = 1;
	params(:,6) = [1 1 -2 2 3 4 5 6 1 1]';
	arg.replace = 0;
end

params(:,6) = params(:,6) * arg.hu_scale;

if streq(arg.type, 'fast') && arg.oversample == 1
	warn('ignoring ''fast'' option for oversample=1')
	arg.type = '';
end

switch arg.type
case 'fast'
	fun = @ellipse_im_fast;
case {'', 'slow'}
	fun = @ellipse_im_slow;
otherwise
	fail('unknown type %s', arg.fast)
end

[phantom, params] = fun(ig.nx, ig.ny, params, ...
	ig.dx, ig.dy, ig.offset_x, ig.offset_y, ...
	arg.rot, arg.oversample, arg.replace);


 function [phantom, params] = ellipse_im_old(nx, ny, params, varargin)
%function [phantom, params] = ellipse_im_old(nx, ny, params, options)
%
% generate ellipse phantom image from parameters:
%	[x_center y_center x_radius y_radius angle_degrees amplitude]
%
% in
%	nx,ny			image size
%	params			ellipse parameters, if empty, use shepp-logan
%
% options:
%	'dx'	?		pixel size
%	'dy'	?		"" (default: -dx for aspire consistency)
%	'fov'	?		use for scaling shepp-logan units
%	'rot'	?		rotate image by this amount [degrees]
%	'oversample'	int	oversampling factor for grayscale boundaries
%	hu_scale	?	use 1000 to scale shepp-logan to HU (default: 1)
%
% note: op ellipse in aspire with nint=3 is oversample=4 = 2^(3-1) here

if ~isvar('ny') || isempty(ny), ny = nx; end

arg.dx = 1;
arg.dy = [];
arg.fov = nx;
arg.rot = 0;
arg.oversample = 1;
arg.hu_scale = 1;
arg = vararg_pair(arg, varargin);

if isempty(arg.dy)
	arg.dy = -arg.dx; % trick: default to match aspire
end

if ~isvar('params') || isempty(params)
	params = shepp_logan_parameters(arg.fov, arg.fov);
end
params(:,6) = params(:,6) * arg.hu_scale;

[phantom, params] = ellipse_im_slow(nx, ny, params, ...
	arg.dx, arg.dy, 0, 0, arg.rot, arg.oversample, 0);


%
% ellipse_im_slow()
% brute force fine grid - can use lots of memory
%
function [phantom, params] = ellipse_im_slow(nx, ny, params, dx, dy, ...
	offset_x, offset_y, rot, over, replace)

if size(params,2) ~= 6
	fail('bad ellipse parameter vector size')
end

% optional rotation of ellipse parameters
if rot ~= 0
	th = deg2rad(rot);
	cx = params(:,1);
	cy = params(:,2);
	params(:,1) = cx * cos(th) + cy * sin(th);
	params(:,2) = -cx * sin(th) + cy * cos(th);
	params(:,5) = params(:,5) + rot;
	clear cx cy th
end

wx = (nx*over-1)/2 + offset_x * over;
wy = (ny*over-1)/2 + offset_y * over;
xx = ((0:nx*over-1) - wx) / over * dx;
yy = ((0:ny*over-1) - wy) / over * dy;
[xx yy] = ndgrid(xx, yy); % fine grid, equally spaced

phantom = zeros(nx*over, ny*over, 'single'); % fine array

ticker reset
ne = nrow(params);
for ie = 1:ne
	ticker(mfilename, ie, ne)

	ell = params(ie, :);
	cx = ell(1);	rx = ell(3);
	cy = ell(2);	ry = ell(4);
	theta = deg2rad(ell(5));
	[xr yr] = rot2(xx-cx, yy-cy, theta);
	tmp = (xr / rx).^2 + (yr / ry).^2 <= 1;

	if replace
		phantom(tmp > 0) = ell(6);
	else
		phantom = phantom + ell(6) * tmp;
	end
end

phantom = downsample2(phantom, over);


%
% ellipse_im_fast()
%
function [phantom, params] = ellipse_im_fast(nx, ny, params, dx, dy, ...
	offset_x, offset_y, rot, over, replace)

if size(params,2) ~= 6
	fail('bad ellipse parameter vector size')
end

% optional rotation
if rot ~= 0
	th = deg2rad(rot);
	cx = params(:,1);
	cy = params(:,2);
	params(:,1) = cx * cos(th) + cy * sin(th);
	params(:,2) = -cx * sin(th) + cy * cos(th);
	params(:,5) = params(:,5) + rot;
	clear cx cy th
end

phantom = zeros(nx, ny, 'single');

wx = (nx-1)/2 + offset_x;
wy = (ny-1)/2 + offset_y;
x1 = ((0:nx-1) - wx) * dx;
y1 = ((0:ny-1) - wy) * dy;
[xx yy] = ndgrid(x1, y1);

if over > 1
	tmp = ((1:over) - (over+1)/2) / over;
	[xf yf] = ndgrid(tmp*dx, tmp*dy);
	xf = xf(:)';
	yf = yf(:)';

	hx = abs(dx) / 2;
	hy = abs(dy) / 2;
end

ticker reset
ne = nrow(params);
for ie = 1:ne
	ticker(mfilename, ie, ne)

	ell = params(ie, :);
	cx = ell(1);	rx = ell(3);
	cy = ell(2);	ry = ell(4);
	theta = deg2rad(ell(5));

	xs = xx - cx; % shift per ellipse center
	ys = yy - cy;

	% coordinates of "outer" corner of each pixel, relative to ellipse center
	xo = xs + sign(xs) * hx; % todo: undefined if over=1
	yo = ys + sign(ys) * hy;

	% voxels that are entirely inside the ellipse:
	[xr, yr] = rot2(xo, yo, theta);
	vi = (xr / rx).^2 + (yr / ry).^2 <= 1;
	gray = single(vi);

	if 0 % grid plotting utility
		pgrid1 = @(x,y,t) plot( ...
			repmat([x(1) x(end)], [length(y) 1])', repmat(y(:), [1 2])', t, ...
			repmat(x(:), [1 2])', repmat([y(1) y(end)], [length(x) 1])', t);
		pgrid = @(x,y,t) pgrid1(x1+hx, y1+hy, t);
	end

	if 0 % examine "inside" voxels
		clf, plot(xx(vi), yy(vi), 'g+')
		hold on, pgrid(x1, y1, 'y-'), hold off
		axis equal, axis square
		plot_ellipse(cx, cy, rx, ry, theta, 'hold', 1)
	end

	if over > 1

		% coordinates of "inner" corner of each pixel, relative to ellipse center
		xi = xs - sign(xs) * hx;
		yi = ys - sign(ys) * hy;

		% voxels that are entirely outside the ellipse:
		[xr yr] = rot2(xi, yi, theta);
		vo = (max(abs(xr),0) / rx).^2 + (max(abs(yr),0) / ry).^2 >= 1;

		if 0 % examine "outside" voxels
			clf, plot(xx(vo), yy(vo), 'ro')
			hold on, pgrid(x1, y1, 'y-'), hold off
			axis equal, axis square
			plot_ellipse(cx, cy, rx, ry, theta, 'hold', 1)
		end

		% subsampling for edge voxels
		edge = ~vi & ~vo;
		x = xx(edge);
		y = yy(edge);
		x = outer_sum(x, xf);
		y = outer_sum(y, yf);

		[xr yr] = rot2(x - cx, y - cy, theta);
		in = (xr / rx).^2 + (yr / ry).^2 <= 1;
		tmp = mean(in, 2);

		if 0 % examine sampling
			clf, plot( ...
				xx(vo), yy(vo), 'ro', ...
				xx(vi), yy(vi), 'g+', ...
				x(in), y(in), 'g.', ...
				x(~in), y(~in), 'r.')
			hold on, pgrid(x1, y1, 'y-'), hold off
			axis equal, axis square
			plot_ellipse(cx, cy, rx, ry, theta, 'hold', 1)
		end

		gray(edge) = tmp;
	end

	if replace
		phantom(gray > 0) = ell(6);
	else
		phantom = phantom + ell(6) * gray;
	end
end


%
% rot2()
% 2d rotation
%
function [xr, yr] = rot2(x, y, theta)
xr = cos(theta) * x + sin(theta) * y;
yr = -sin(theta) * x + cos(theta) * y;


%
% parameters from Kak and Slaney text, p. 255
% the first four columns are unitless "fractions of field of view"
%
function params = shepp_logan_parameters(xfov, yfov)
params = [...
	0	0	0.92	0.69	90	2;
	0	-0.0184	0.874	0.6624	90	-0.98;
	0.22	0	0.31	0.11	72	-0.02;
	-0.22	0	0.41	0.16	108	-0.02;
	0	0.35	0.25	0.21	90	0.01;
	0	0.1	0.046	0.046	0	0.01;
	0	-0.1	0.046	0.046	0	0.01;
	-0.08	-0.605	0.046	0.023	0	0.01;
	0	-0.605	0.023	0.023	0	0.01;
	0.06	-0.605	0.046	0.023	90	0.01];

params(:,[1 3]) = params(:,[1 3]) * xfov/2;
params(:,[2 4]) = params(:,[2 4]) * yfov/2;


%
% ellipse_im_profile()
%
function ellipse_im_profile
ig = image_geom('nx', 2^9, 'ny', 2^9+2', 'fov', 250);
profile on
x0 = ellipse_im(ig, [], 'oversample', 3, 'type', 'fast');
profile off
profile report


% ellipse_im_test()
function ellipse_im_test
ig = image_geom('nx', 2^8, 'ny', 2^8+2', 'fov', 250);
%ig.offset_y = 75.6 / ig.dy;

im plc 2 3

over = 2^2;
x0 = ellipse_im(ig, [], 'oversample', over, 'fov', 250);
im(1, ig.x, ig.y, x0, 'Shepp Logan', [0.9 1.1]), cbar

x1 = ellipse_im(ig, 'shepplogan-emis', 'oversample', over, 'fov', 250);
im(4, x1, 'Shepp Logan Emission'), cbar

x3 = ellipse_im(ig, 'shepplogan-brainweb', 'fov', 250);
im(5, x3, 'Shepp Logan Brainweb'), cbar


if 0 % test vs old shepplogan
	x2 = shepplogan(ig.nx, ig.ny, 1);
	im(3, x1), cbar
	im(4, x2-x1), cbar
return
end

% compare to aspire
if ~has_aspire, return, end

nx = 2^6;
ig = image_geom('nx', nx, 'ny', nx+2, 'fov', 2^7);
ell = [10 20 30 40 50 1];

file = [test_dir '/t.fld'];
com = sprintf('echo y | op -chat 99 ellipse %s %d %d  %g %g %g %g %g %g %d', ...
	file, ig.nx, ig.ny, ell ./ [ig.dx ig.dx ig.dx ig.dx 1 1], log2(over)+1);
os_run(com)
asp = fld_read(file);
im(4, asp, 'aspire'), cbar

area.asp = sum(asp(:)) * abs(ig.dx * ig.dy);

types = {'slow', 'fast'};
for ii=1:length(types)
	type = types{ii};
	cpu etic
	mat = ellipse_im(ig, ell, 'oversample', over, 'type', types{ii});
	area.(type) = sum(mat(:)) * abs(ig.dx * ig.dy);
	cpu('etoc', [type ' time:'])

	im(1+ii, mat, sprintf('mat, %s', type)), cbar
	im(4+ii, (mat-asp)*over^2, 'difference: (mat-asp)*over$^2$'), cbar

	pr max(abs(col(mat - asp))) / ell(6) * over^2

%	equivs(mat, asp)
%	max_percent_diff(mat, asp)
end

area.real = pi * ell(3) * ell(4) * ell(6);
pr area
