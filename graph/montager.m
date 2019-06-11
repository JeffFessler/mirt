 function xo = montager(xi, varargin)
%function xo = montager(xi, varargin)
%|
%| Make montage or mosaic of set of images.
%|
%| in
%|	xi		[3d or 4d] set of images
%|
%| options
%|	'col'		# of cols
%|	'row'		# of rows
%|	'aspect'	aspect ratio (default 1.2)
%|	'byte'		if 1, convert each slice to [0,255]
%|
%| out
%|	xo [2d]		3d or 4d images arrange a as a 2d rectangular montage
%|
%| Copyright 1997, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(xi, 'test'), montager_test, return, end

arg.aspect = 1.2; % trick: default 1.2/1 aspect ratio
arg.col = [];
arg.row = [];
arg.byte = false;
arg.chat = 0;

arg = vararg_pair(arg, varargin);

dims = size(xi);
nx = dims(1);
ny = dims(2);
nz = prod(dims(3:end));
if ndims(xi) > 3
	xi = reshape(xi, [nx ny nz]);
end

if isempty(arg.col)
	if isempty(arg.row)
		if ndims(xi) == 4
			arg.col = n3;
		elseif nx == ny && nz == round(sqrt(nz)).^2 % perfect square
			arg.col = round(sqrt(nz));
		else
			arg.col = ceil(sqrt(nz * ny / nx * arg.aspect));
		end
	else
		arg.col = ceil(nz / arg.row);
	end
end

if isempty(arg.row)
	arg.row = ceil(nz / arg.col);
end

if arg.chat
	pr '[arg.row arg.col]'
end

if arg.byte
	xo = zeros(nx * arg.col, ny * arg.row, 'uint8');
else
	xo = zeros(nx * arg.col, ny * arg.row);
end

for iz=0:(nz-1)
	iy = floor(iz / arg.col);
	ix = iz - iy * arg.col;
	tmp = xi(:,:,iz+1);
	if arg.byte
		scale = div0(255, max(tmp(:)) - min(tmp(:)));
		tmp = floor((tmp-min(tmp(:))) * scale);
	end
	xo([1:nx]+ix*nx, [1:ny]+iy*ny) = tmp;
end


% montager_test()
function montager_test
t = [20 30 5];
%t = [224 4422 8];
t = reshape([1:prod(t)], t);
im plc 1 3
im(1, montager(t, 'chat', 0))
im(2, montager(t, 'row', 4))
im(3, montager(t, 'byte', 1, 'row', 1))
