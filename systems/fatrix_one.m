  function ob = fatrix_one(dim, fun, arg, varargin)
%|function ob = fatrix_one(dim, fun, arg, varargin)
%|
%| Construct Fatrix object where one column of it is computed at a time by 'fun'
%|
%| in
%|	dim	[1 2]	Fatrix dimension
%|	fun	@	fun(arg, n) must return nth column of matrix
%|	arg	struct	argument for fun
%|
%| option
%|	'type'	char	'col' only (default)
%|	'par'	0|1	if nonzero, use 'parfor' (default: true)
%|			caution: using par=1 can make the results of A*x
%|			vary slightly between runs because the accumulation
%|			operation depends very slightly on execution order
%|
%| out
%|	ob		Fatrix
%|
%| Copyright 2007-12-14, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(dim, 'test'), fatrix_one_test, return, end
if nargin == 1 && streq(dim, 'test0'), fatrix_one_test0, return, end

arg.dim = dim;
arg.fun = fun;
arg.arg = arg;

% options
arg.type = 'col';
arg.par = true;
arg.precise = 'single';
arg = vararg_pair(arg, varargin);

if length(arg.dim) ~= 2, fail 'bad dim', end

% build Fatrix object
switch arg.type
case 'col'
	ob = Fatrix(arg.dim, arg, 'caller', 'fatrix_one(col)', ...
		'forw', @fatrix_one_col_forw, 'back', @fatrix_one_col_back);
otherwise
	fail 'not done'
end


%
% fatrix_one_col_forw(): y = G * x
%
function y = fatrix_one_col_forw(arg, x)

dimx = size(x);
diml = dimx(2:end);
x = reshape(x, arg.dim(2), prod(diml)); % [N *L]

y = zeros(arg.dim(1), prod(diml), arg.precise); % [M *L]
if arg.par
	parfor nn=1:arg.dim(2)
		tmp = arg.fun(arg, nn);
		y = y + tmp * x(nn,:);
	end
else
	for nn=1:arg.dim(2)
		tmp = arg.fun(arg, nn);
		y = y + tmp * x(nn,:);
	end
end
y = reshape(y, arg.dim(1), diml); % [M (L)]


%
% fatrix_one_col_back(): x = G' * y
%
function x = fatrix_one_col_back(arg, y)

dimy = size(y);
diml = dimy(2:end);
y = reshape(y, arg.dim(1), prod(diml)); % [M *L]

x = zeros(arg.dim(2), prod(diml), arg.precise); % [N *L]
if arg.par
	parfor nn=1:arg.dim(2)
		tmp = arg.fun(arg, nn);
		x(nn,:) = tmp' * y;
	end
else
	for nn=1:arg.dim(2)
		tmp = arg.fun(arg, nn);
		x(nn,:) = tmp' * y;
	end
end
x = reshape(x, arg.dim(2), diml); % [N (L)]


% fatrix_one_tester()
function fatrix_one_tester(par)

dim = [4 5];
A1 = reshape(prod(dim):-1:1, dim);
tmp.store = A1;
fun = @(tmp, nn) tmp.store(:,nn);

A2 = fatrix_one(dim, fun, tmp, 'par', par);

x = [1:dim(2)]';
y1 = A1 * x;
y2 = A2 * x;
if par
	equivs(y1, y2)
else
	jf_equal(y1, y2)
end

x1 = A1' * y1;
x2 = A2' * y1;
if par
	equivs(y1, y2)
else
	jf_equal(x1, x2)
end

Fatrix_test_basic(A2, true(dim(2),1))


% fatrix_one_test0()
% test without parallel toolbox
function fatrix_one_test0
fatrix_one_tester(0)


% fatrix_one_test()
function fatrix_one_test

for par=0:1
	fatrix_one_tester(par)
end
