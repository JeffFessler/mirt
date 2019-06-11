 function y = reshapee(x, varargin)
%|function y = reshapee(x, varargin)
%|
%| reshape function that allows possibly one null argument, and all
%| other arguments can be vectors, unlike matlab that requires scalars.
%| example: reshape(rand(2*3*5,7), [2 3], [], 7) will become [2 3 5 7]
%|
%| in
%|	x		[(*dim)]
%|	varargin	dimensions
%|
%| out
%|	y		[dim]
%|
%| Copyright 2004-8-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), reshapee_test, return, end
if nargin < 2, ir_usage, end

dim_i = size(x);

ndim = 0;
edim = [];
dim_o = [];
for ii=1:length(varargin)
	arg = varargin{ii};
	ndim = ndim + length(arg);
	if isempty(arg)
		if ~isempty(edim)
			fail 'only one empty dim allowed'
		end
		edim = 1+length(dim_o);
		dim_o = [dim_o, 1]; % trick: place holder
	else
		dim_o = [dim_o, arg];
	end
end

if ~isempty(edim) % fill in empty dim if present
	if prod(dim_o) <= 0, fail('nonpositive dims?'), end
	dim_o(edim) = prod(dim_i) / prod(dim_o);
	if round(dim_o(edim)) ~= dim_o(edim)
		pr dim_i
		pr dim_o
		fail('bad dim')
	end
end

y = reshape(x, dim_o);


% self test
function reshapee_test
dim = [2 3 5 7];
x = 1:prod(dim);
y = reshapee(x, dim(1:2), [], dim(4));
z = reshape(x, dim);
jf_equal(y,z)
