 function y = ir_dct2(x)
%function y = ir_dct2(x)
%|
%| 2D DCT via separable 1D processing for many 2D images
%| because matlab dct2() only handles a single 2D image.
%| in
%| x [nx ny npatch]
%| out
%| y [nx ny npatch]
%| 
%| 2017 Jeff Fessler, University of Michigan

if nargin ~= 1, ir_usage, end

[nx, ny, np] = size(x);
y = permute(reshape(...
	dct(reshape(permute(reshape(...
	dct(reshape(x, nx, [])), nx, ny, np), ...
	[2 1 3]), ny, [])), ny, nx, np), [2 1 3]);
