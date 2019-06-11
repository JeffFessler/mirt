  function ss = outer_sum(xx,yy)
%|function ss = outer_sum(xx,yy)
%|
%| compute an "outer sum" x + y'
%| that is analogous to the "outer product" x * y'
%|
%| in
%|	xx	[nx 1]
%|	yy	[1 ny]
%|		more generally: xx [(dim)] + yy [L,1] -> xx [(dim) LL]
%| out
%|	ss [nx ny]	ss(i,j) = xx(i) + yy(j)
%|
%| Copyright 2001, Jeff Fessler, University of Michigan

if ~nargin, ir_usage, end
if streq(xx, 'test'), outer_sum_test, return, end

% for 1D vectors, allow rows or cols for backward compatibility
if ndims(xx) == 2 && min(size(xx)) == 1 && ndims(yy) == 2 && min(size(yy)) == 1
	nx = length(xx);
	ny = length(yy);
	xx = repmat(xx(:), [1 ny]);
	yy = repmat(yy(:)', [nx 1]);
	ss = xx + yy;
return
end

%if size(xx,1) == 1
%	warn 'xx is a row vector? are you sure?'
%end

% otherwise, xx is not a vector, but yy must be

xdim = size(xx);
ydim = size(yy);
if ndims(yy) ~= 2 || min(size(yy)) ~= 1
	error 'yy must be a vector'
end
sdim = [xdim length(yy)];
xo = repmat(xx, [ones(1, length(xdim)) length(yy)]); % [xdim] -> [xdim ny]
yo = repmat(yy(:), [1 xdim]); yo = permute(yo, [2 3 1]);
% yo = repmat(yy(:)', [xdim 1 1]); % not sure how to make this work.
ss = xo + yo;

function outer_sum_test
%xx = [1:4];
%yy = [0:10:50];
%pr outer_sum(xx,yy)
pr 'outer_sum([1:4], [0:10:50])'
xx = outer_sum([1:4], [0:10:50]);
pr 'outer_sum(xx, [100 200])'
