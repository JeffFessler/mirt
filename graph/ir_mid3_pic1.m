  function b = ir_mid3_pic1(x, file, varargin)
%|function b = ir_mid3_pic1(x, file, varargin)
%|
%| write 'mid3' type picture to .png file
%|
%| in
%|	x	[nx ny nz]	3D image
%|	file	char		file name
%|
%| option
%|	'clim'			display limit, default [800 1200]
%|
%| out
%|	b	[nx+nz ny+nz]	2D 'mid3' image, uint8
%|
%| 2012-10-01, Jeff Fessler

if nargin < 1, help(mfilename), error(mfilename), end
if streq(x, 'test'), ir_mid3_pic1_test, return, end

arg.clim = [800 1200];
arg = vararg_pair(arg, varargin);

y = jf_mip3(x, 'type', 'mid');

clim = arg.clim;
b = 255 * (y - clim(1)) / (clim(2) - clim(1));
b = min(b,255);
b = max(b,0);
b = uint8(b);

if exist(file, 'file')
	warn('file "%s" exists already.', file)
	yn = input('over-write? y|n ', 's');
	if ~streq(yn, 'y', 1)
		error('file name')
	end
end
imwrite(b', file); % trick: transpose


% ir_mid3_pic1_test
function ir_mid3_pic1_test

ig = image_geom('nx', 64, 'nz', 30, 'dx', 1, 'dz', 2);
x = ellipsoid_im(ig, [], 'hu_scale', 1000);

file = 'test.png';
b = ir_mid3_pic1(x, file, 'clim', [950 1050]);

%h = im(y, arg.clim);
%z = get(h, 'cdata');
%minmax(b)

tmp = imread(file)'; % trick
jf_equal(tmp, b)

% !xv test.png &
