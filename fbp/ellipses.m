 function [phantom, params] = ellipses(nx, ny, params, dx, dy, varargin)
%function [phantom, params] = ellipses(nx, ny, params, dx, dy, options)
%
% generate ellipse phantom image from parameters:
% [x_center y_center x_radius y_radius angle_degrees amplitude [oversample]]
%
% in
%	nx,ny			image size
%	params			ellipse parameters, if empty, use shepp-logan
%	(dx,dy)			pixel size
%
% options:
%	'fov' (value)		use for scaling shepp-logan units
%	'erot' (value)		rotate all ellipses by this amount [degrees]
%	'oversample'		specify oversampling factor here instead
%
% outdated options:
%	'xscale' (value)	use -1 to flip in x (or use negative dx)
%	'yscale' (value)	use -1 to flip in y (or use negative dy)
%
% note: op ellipse in aspire with nint=3 is oversample=4 = 2^(3-1) here
%
% Copyright 2003-5-16, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

if nargin == 1 && streq(nx, 'test'), ellipses_test, return, end

if nargin == 1 && streq(nx, 'shepplogan128')
	phantom = ellipses(2^7);
return
end

warning 'ellipses() has been replaced by ellipse_im()'

if ~isvar('ny') || isempty(ny)
	ny = nx;
end

if ~isvar('dx') || isempty(dx)
	dx = 1;
end

if ~isvar('dy') || isempty(dy)
	dy = dx;
end

arg.xscale = 1;
arg.yscale = 1;
arg.fov = nx;
arg.erot = 0;
arg.oversample = 1;
arg = vararg_pair(arg, varargin);

if ~isvar('params') || isempty(params)
	params = shepp_logan_parameters(arg.fov, arg.fov);
end

% optional rotation
if arg.erot ~= 0
	th = deg2rad(arg.erot);
	x = params(:,1);
	y = params(:,2);
	params(:,1) = x * cos(th) + y * sin(th);
	params(:,2) = -x * sin(th) + y * cos(th);
	params(:,5) = params(:,5) + arg.erot;
	clear x y
end

phantom = zeros(nx, ny);

ticker reset
ne = nrow(params);
for ie = 1:ne;
	ticker(mfilename, ie, ne)

	ell = params(ie, :);

	if length(ell) == 6
		over = arg.oversample;
	elseif length(ell) == 7
		over = ell(7);
		warning 'oversampling in ell is deprecated.  use option instead'
	else
		error 'bad ellipse parameter vector size'
	end

	xx = arg.xscale * ([1:nx*over] - (nx*over+1)/2) * dx / over;
	yy = arg.yscale * ([1:ny*over] - (ny*over+1)/2) * dy / over;
	[xx yy] = ndgrid(xx, yy);

	cx = ell(1);		rx = ell(3);
	cy = -sign(dy)*ell(2);	ry = ell(4);
	theta = -ell(5) / 180 * pi;% * sign(arg.yscale);
	x = cos(theta) * (xx-cx) + sin(theta) * (yy-cy);
	y = -sin(theta) * (xx-cx) + cos(theta) * (yy-cy);
	tmp = (x/rx).^2 + (y/ry).^2 <= 1;

	phantom = phantom + ell(6) * downsample2(tmp, over);
end


%
% parameters from Kak and Slaney text, p. 255
% most of these values are unitless "fractions of field of view"
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
function ellipses_test
x = ellipses(2^7);
im(x, 'Shepp Logan', [0.9 1.1]), cbar

if ~has_aspire, return, end

% compare to aspire
nx = 2^7; ny = nx+2;
ell = [10 20 30 40 50 100];
mat = ellipses(nx, ny, ell, 1, 1, 'oversample', 4);
dir = test_dir;
file = [dir '/t.fld'];
com = sprintf('echo y | op ellipse %s %d %d  %g %g %g %g %g %g 3', ...
	file, nx, ny, ell)
os_run(com)
asp = fld_read(file);

im clf, pl=220;
im(pl+1, mat, 'mat'), cbar
im(pl+2, asp, 'aspire'), cbar
im(pl+3, asp-mat, 'difference'), cbar
max_percent_diff(mat, asp)
