 function y = ir_phantom_resize(x, nx, ny)
%function y = ir_phantom_resize(x, nx, ny)
%| resize a [mx,my] phantom image to be [nx,ny]
%| by combination of downsampling and cropping

if nargin == 1 && streq(x, 'test'), ir_phantom_resize_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

mx = size(x,1);
my = size(x,2);
if nx == mx && ny == my
	y = x;
	return
end

if nx > mx || ny > my
	error 'requested size too large, implement interpolation!'
end

% if power of 2 size, then downsample
if 2^floor(log2(nx)) == nx && nx == ny
	y = downsample2(x, mx/nx);
return
end

nn = 2^max(ceil(log2([nx ny])));
y = downsample2(x, mx/nn);

% trim if non power of 2
if nx < nn
	y = y([1:nx]+round((nn-nx)/2),:);
end
if ny < nn
	y = y(:,[1:ny]+round((nn-ny)/2));
end


function ir_phantom_resize_test
x = zeros(64, 60);
x(20:30,10:15) = 1;
y = ir_phantom_resize(x, 32, 28);
im plc 1 2
im(1, x)
im(2, y)
