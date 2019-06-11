 function ob = Gpermute(Nd, varargin)
%function ob = Gpermute(Nd, varargin)
%|
%| Construct Gpermute object, which permutes input image dimensions.
%| This is useful perhaps for converting xyz to zyx order, for example.
%| See Gpermute_test() for example usage.
%|
%| in
%|	Nd	[1 D]		input signal dimensions
%|
%| options
%|	mask	logical		default: [] which means true(Nd)
%|	how	[1 D]		e.g. [3 1 2] for xyz to zxy (default: 1:D)
%|
%| out
%|	ob			fatrix2 object
%|
%| Copyright 2012-06-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(Nd, 'test'), Gpermute_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% defaults
arg.mask = [];
arg.idim = Nd;
arg.how = 1:numel(Nd);

% options
arg = vararg_pair(arg, varargin);

jf_equal(sort(arg.how), 1:numel(Nd))
arg.odim = arg.idim(arg.how);
if isempty(arg.mask)
	arg.omask = [];
else
	arg.omask = permute(arg.mask, arg.how);
end

%forw = @(arg, x) permute(x, arg.how);
%back = @(arg, y) ipermute(y, arg.how);
forw = @(arg, x) fatrix2_maskit(arg.omask, permute(x, arg.how));
back = @(arg, y) fatrix2_maskit(arg.mask, ipermute(y, arg.how));
ob = fatrix2('arg', arg, 'omask', arg.omask, 'odim', arg.odim, ...
		'abs', @(ob) ob, 'forw', forw, 'back', back);


% Gpermute_test
function Gpermute_test

dim = [4 5 6];
how = [3 1 2];
mask = true(dim);
mask(end-2) = false;
A = Gpermute(dim, 'how', how, 'mask', mask);

x = reshape(1:prod(dim), dim);

y1 = A.arg.omask .* permute(x, how);
y2 = A * x;
jf_equal(y1, y2)

y2 = A * x(mask(:));
y2 = embed(y2, A.arg.omask);
jf_equal(y1, y2)

y = y1;
x1 = ipermute(y, how);
x2 = A' * y;
jf_equal(x1, x2)

x2 = A' * y(A.arg.omask);
x2 = embed(x2, mask);
jf_equal(x1, x2)

fatrix2_tests(A)
test_adjoint(A);
