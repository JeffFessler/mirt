 function [out weight] = ir_patch_avg(patches, idx, dim, varargin)
%function [out weight] = ir_patch_avg(patches, idx, dim, varargin)
%|
%| average 2D overlapping patches to form image
%| each output pixel value is average of corresponding pixel values
%| from all patches that contain that pixel, i.e, that overlapping it
%|
%| in
%|	patches	[*patch_size npatch]
%|	idx	[npatch]		from ir_im2col
%|
%| option
%|	patch_size	[2]		patch size (default: [8 8])
%|	means	[npatch]		patch means (if need subtracted)
%|						default: empty
%|	nchunk	scalar			# of patches to due in block processing
%|						default: 10000
%|
%| out
%|	out	[dim]	 image formed by averaging overlapping patches
%|
%| based on code from Sai Ravishankar
%| 2016-03-03, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(patches, 'test'), ir_patch_avg_test, return, end

arg.patch_size = [8 8];
arg.means = [];
arg.nchunk = 10000; % how many blocks of patches to do jointly
arg = vararg_pair(arg, varargin);

npatch = size(patches, 2);

if numel(arg.patch_size) ~= 2 || (prod(arg.patch_size) ~= size(patches,1))
	fail 'bad patch_size'
end

b1 = arg.patch_size(1);
b2 = arg.patch_size(2);

[rows, cols] = ind2sub(dim - [b1 b2] + 1, idx);

if ~isempty(arg.means)
	if numel(arg.means) ~= npatch, fail 'bad arg.means', end
	patches = patches + repmat(arg.means, [b1*b2 1]);
end

out = zeros(dim);
weight = zeros(dim);

for jj = 1:arg.nchunk:npatch
	jumpSize = min(jj+arg.nchunk-1, npatch);
	zz = patches(:, jj:jumpSize);
	for ii = jj:jumpSize
		col = cols(ii); row = rows(ii);
		block = reshape(zz(:, ii-jj+1), [b1 b2]);
		i1 = row:(row+b1-1);
		i2 = col:(col+b2-1);
		out(i1, i2) = out(i1, i2) + block; % +=
		weight(i1, i2) = weight(i1, i2) + 1;
	end
end

if any(weight(:) == 0)
	fail 'bug'
end
out = out ./ weight; % average


% ir_patch_avg_test
function ir_patch_avg_test

nx = 2^8 - 8*7*0; ny = 2^8;
xtrue = shepplogan(nx, ny, 1);
%xtrue = rand(nx, ny);
xtrue(end) = 5; % stress ends
xtrue(nx) = 4; % stress ends
bsize = [8 6];
stride = 3;
stride = 1;
[patches idx] = ir_im2col(xtrue, bsize, stride);
im plc 2 2
im(1, xtrue)
[xhat weight] = ir_patch_avg(patches, idx, [nx ny], 'patch_size', bsize);
im(2, xhat)
im(3, weight)
im(4, xhat - xtrue)
equivs(xhat, xtrue)
