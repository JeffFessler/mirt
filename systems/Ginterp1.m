 function ob = Ginterp1(mask, varargin)
%function ob = Ginterp1(mask, options)
%|
%| Construct Ginterp1 object for 1D image registration.
%| This method does not have an adjoint.
%|
%| See Ginterp1_test() below for example usage.
%|
%| in
%|	mask	size(image)	logical array of object support.
%|
%| options
%|	'interp1_arg'	{}	arguments to interp1 for 'forward' operation
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|			nd = numel(interp1_arg{2})
%|
%| Copyright 2013-12-05, Jeff Fessler, University of Michigan

if nargin == 1 && streq(mask, 'test'), Ginterp1_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = mask;

% option defaults
arg.interp1_x = [];
arg.interp1_xq = [];
arg.interp1_arg = {'linear', 0};

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

if isempty(arg.interp1_x) || isempty(arg.interp1_xq)
	fail '"interp1_arg" and "interp1_x" required'
end

arg.ndim = ndims(mask);
if arg.ndim == 2 && size(mask,2) == 1
	arg.ndim = 1;
else
	fail 'only 1D mask allowed'
end

if streq(arg.interp1_arg{1}, 'linear')
	abs_arg = {'abs', @(ob) ob}; % linear interpolation uses nonnegative coefficients
else
	abs_arg = {}; % unknown
end

arg.fun_forw = @(arg, x) ...
	interp1(arg.interp1_x, x, arg.interp1_xq, arg.interp1_arg{:});
arg.fun_back = @(arg, y) fail('interp adjoint not done');

% build object
idim = size(mask);
if numel(idim) == 2 && idim(2) == 1
	idim = idim(1); % 1d
end
ob = fatrix2('mask', mask, 'arg', arg, ...
	'idim', idim, 'odim', numel(arg.interp1_xq), ...
	abs_arg{:}, 'forw', arg.fun_forw, 'back', arg.fun_back);


% Ginterp1_test()
function Ginterp1_test

nx = 8;
mask = true(nx,1);
mask(1) = false;
xq = linspace(0,nx+1, 101);

A = Ginterp1(mask, 'interp1_x', 1:nx, 'interp1_xq', xq);

if 0 % todo
	fatrix2_tests(A, 'complex', 0, 'halt', 0, ...
		'check1', false, 'full', false) % because of bad adjoint
%	test_adjoint(A, 'complex', 1);
else
	warn 'adjoint not tested'
end

im plc 1 2
Af = full(A);
im(1, Af'), axis xy, axis normal

im subplot 2
x = [1:nx]' .* mask;
y = A * x;

plot(1:nx, x, 'o', xq, y, '.-')
