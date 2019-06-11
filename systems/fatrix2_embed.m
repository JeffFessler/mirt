 function z = fatrix2_embed(mask, dim, x)
%function z = fatrix2_embed(mask, dim, x)
%|
%| version of embed that handles an empty mask
%|
%| in
%|	x	[N L]
%|	dim	[]		size(mask) if mask is nonempty
%|	mask	[(dim)]		possibly empty
%|
%| out
%|	z	[(dim) L]
%|
%| Copyright 2010-12-21, Jeff Fessler, University of Michigan

if nargin < 3, ir_usage, end

if ~isempty(mask)
	z = embed(x, mask);
else
	z = reshape(x, [dim size(x,2)]);
end
