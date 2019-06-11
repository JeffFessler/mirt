  function y = reshaper(x, dim)
%|function y = reshaper(x, dim)
%|
%| reshape function that is more flexible, allowing for "multiples".
%| example: reshape(rand(2*3*5,7), [2 3 5]) will become [2 3 5 7]
%|
%| in
%|	x	[*dim (Ld)]
%|	dim	short row or column
%|		if dim is '2d' then y is 2d with first n-1 dims collapsed
%|
%| out
%|	y	[dim (Ld)]
%|
%| Copyright 2004-8-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), reshaper_test, return, end
if nargin < 2, ir_usage, end

dim_i = size(x);

if ischar(dim) && streq(dim, '2d')
	y = reshape(x, prod(dim_i(1:end-1)), dim_i(end));
return
end

dim_e = dim_i(2:end); % extra dimensions (this could be made fancier)
dim_o = [dim dim_e];
y = reshape(x, dim_o);


% self test
function reshaper_test
x = reshape(1:2*3*5*7, 2*3*5, 7);
dim = [2 3 5];
y = reshaper(x, dim);
jf_equal(size(y), [dim 7])
