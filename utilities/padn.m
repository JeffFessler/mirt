  function out = padn(mat, newdim)
%|function out = padn(mat, newdim)
%| pad input matrix to newdim, preserving 'center point'
%| using this avoids matlab's padarray.m which is in image proc toolbox

if nargin < 1, ir_usage, end
if streq(mat, 'test'), padn_test, return, end

% default: round to next up power of 2
if nargin < 2
	newdim = 2 .^ ceil(log2(size(mat)))
end

olddim = size(mat);
if any(newdim < olddim), error('must be bigger'), end

idx = cell(1, length(newdim));
for ii=1:length(newdim)
	pad = newdim(ii) - olddim(ii);
	if ~rem(pad,2) % even
		offset = pad/2;
	else
		if rem(olddim(ii),2) % odd
			offset = ceil(pad/2);
		else
			offset = floor(pad/2);
		end
	end
	idx{ii} = [1:olddim(ii)] + offset;
end
out = zeros1(newdim);
out(idx{:}) = mat;


 function z = zeros1(dim)
%function z = zeros1(dim)
% make array of zeros(dims(1), dims(2), ...)
% works logically even if the input is a scalar
% jeff fessler

if nargin < 1, ir_usage, end

if length(dim) == 1
	dim = [dim 1]
end
z = zeros(dim);


function padn_test
x = ones(3,4);
jf_equal(x, unpadn(padn(x, [4 6]), size(x)))
jf_equal(x, unpadn(padn(x, [5 7]), size(x)))

jf_equal(padn([1 2 1], [1 5]), [0 1 2 1 0])
jf_equal(padn([0 1 2 1], [1 5]), [0 1 2 1 0])
jf_equal(padn([0 1 2 1], [1 6]), [0 0 1 2 1 0])
jf_equal(padn([1 2 1], [1 6]), [0 0 1 2 1 0])
