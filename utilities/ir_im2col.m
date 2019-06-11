 function [blocks, idx] = ir_im2col(I, blkSize, stride)
%function [blocks, idx] = ir_im2col(I, blkSize, stride)
%|
%| Extract patches of 2D image and corresponding index of first "corner" of each
%| Variant of 'im2col' that allows arbitrary stride between patches.
%| Reverts to 'im2col' if stride == 1.
%|
%| in
%|	I	[nx ny]			image from which to extract patches
%|	blkSize	[1]|[2]			patch size
%|	stride	[1]|[2]			sliding distance
%|
%| out
%|	blocks	[*blkSize npatch]	extracted patches
%|					npatch approx nx/stride(1)*ny/stride(2)
%|	idx	[npatch]		index of patch start
%|						caution: tricky indexing
%|
%| Original version from K-SVD code
%| http://www.cs.technion.ac.il/~elad/software/
%| http://www.ifp.illinois.edu/~yoram/DLMRI-Lab/Documentation.html
%|
%| 2016-03-03 Jeff Fessler, University of Michigan, add documentation, etc.

if nargin < 1, ir_usage, end
if nargin == 1 && streq(I, 'test'), ir_im2col_test, return, end

switch numel(blkSize)
case 1
	blkSize = [blkSize blkSize];
case 2
	% ok
otherwise
	fail('bad blkSize')
end

switch numel(stride)
case 1
	stride = [stride stride];
case 2
	% ok
otherwise
	fail('bad stride')
end

if (all(stride == 1))
	blocks = im2col(I, blkSize, 'sliding');
	idx = [1:size(blocks,2)];
return
end

if any(stride == -1) % trick to check this code vs matlab's im2col
	stride = [1 1];
end

% take blocks (patches) in 'stride' distances,
% but always take the first and last one (in each row and column)
% so that every pixel is in at least one patch
s1 = stride(1);
s2 = stride(2);
idxMat = zeros(size(I) - blkSize + 1);
idxMat([[1:s1:end-1],end], [[1:s2:end-1],end]) = 1;

idx = find(idxMat); % caution: tricky since idxMat is smaller than "I"
[rows, cols] = ind2sub(size(idxMat), idx);
blocks = zeros([blkSize numel(idx)]);
b1 = blkSize(1) - 1; % trick
b2 = blkSize(2) - 1;

% todo: The Mathworks "Hankel indexing" in im2col.m is probably more efficient
% than this brute-force loop approach.  But this is simpler...
for ii = 1:numel(idx)
	blocks(:,:,ii) = I(rows(ii):(rows(ii)+b1), cols(ii):(cols(ii)+b2));
end
blocks = reshape(blocks, prod(blkSize), []);


% ir_im2col_test
function ir_im2col_test
n1 = 4;
n2 = 5;
I = reshape(1:(n1*n2), n1, n2);
b1 = 2;
b2 = 3;
bsize = [b1 b2];
[b, i] = ir_im2col(I, bsize, 2);
b = reshape(b, b1, b2, []);

jf_equal(b(:,:,1), I(1:b1, 1:b2))
jf_equal(i(1), 1)

e1 = n1-b1+1;
e2 = n2-b2+1;
jf_equal(b(:,:,end), I(e1:n1, e2:n2))
jf_equal(i(end), e1*e2)

rng(0)
nx = 2^6; ny = nx - 2;
f1 = rand(nx, ny);

if 1 % test stride=1 case
	p1 = ir_im2col(f1, bsize, 1);
	p2 = ir_im2col(f1, bsize, -1);
	jf_equal(p1, p2)
end

if 1 % compare convolution with patch operation as sanity check
	h = eye(5) + ones(5);
	g1 = conv2(f1, h, 'valid');
	bsize = size(h);
	p2 = ir_im2col(f1, bsize, 1);
	tmp = h(:)' * p2; % should be equivalent to convolution
	g2 = reshape(tmp, nx - bsize(1) + 1, ny - bsize(2) + 1);
	equivs(g1, g2)
end
