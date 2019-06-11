  function ob = Godwt1(mask, varargin)
%|function ob = Godwt1(mask, varargin)
%| Construct Godwt1 object that computes orthonormal discrete wavelet
%| decomposition of a signal with dimensions [(N)].
%| This is useful for sparsity regularization (aka compressed sensing).
%| (1D or 2D wavelets only)
%|
%| in
%|	'mask'	logical [(Nd)]	image-domain mask, often true(nx,ny)
%|
%| options
%|	'level'	int	decomposition level (default 1)
%|	'wname'	char	wavelet name. default: 'haar'
%|
%| out
%|	ob	[*N np]	fatrix  object, where np = sum(mask(:))
%|
%|
%| Copyright 2012-05-17, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(mask, 'test'), Godwt1_test, return, end

arg.mask = mask;
arg.level = 1;
arg.wname = 'haar';
arg.abs = false;
arg = vararg_pair(arg, varargin);

if isempty(arg.mask), fail 'must provide a mask', end

% transform dimension
idim = size(arg.mask);
if idim(end) == 1
	idim = idim(1:end-1);
end

switch numel(idim)
case 1
	forw = @(arg, x) ir_odwt1(x, ...
		'level', arg.level, 'wname', arg.wname, 'abs', arg.abs);
	back = @(arg, y) fatrix2_maskit(arg.mask, ir_odwt1(y, 'adj', 1, ...
		'level', arg.level, 'wname', arg.wname, 'abs', arg.abs));
	does_many = true;
case 2
	forw = @(arg, x) ir_odwt2(x, ...
		'level', arg.level, 'wname', arg.wname, 'abs', arg.abs);
	back = @(arg, y) fatrix2_maskit(arg.mask, ir_odwt2(y, 'adj', 1, ...
		'level', arg.level, 'wname', arg.wname, 'abs', arg.abs));
	does_many = false;
otherwise
        fail('dim %d unsupported', numel(idim))
end

% build fatrix2 object
arg.idim = idim;
ob = fatrix2('idim', idim, 'odim', idim, 'arg', arg, ...
	'abs', @Godwt1_abs, 'meth', {'codes', @Godwt1_codes, '()'}, ...
	'forw', forw, 'back', back, 'does_many', does_many);


% Godwt1_abs()
function ob = Godwt1_abs(ob)
ob.arg.abs = true;


% Godwt1_codes()
function codes = Godwt1_codes(arg)

switch numel(arg.idim)
case 1
	[dummy codes] = ir_odwt1(zeros([arg.idim 1]), ...
		'level', arg.level, 'wname', arg.wname);
case 2
	[dummy codes] = ir_odwt2(zeros(arg.idim), ...
		'level', arg.level, 'wname', arg.wname);
otherwise
        fail('dim %d unsupported', numel(idim))
end


% Godwt1_test()
function Godwt1_test

if 1 % 1d
	mask = true(8*3,1);
	mask(1:3) = false;
	level = 3;
	U = Godwt1(mask, 'level', level);
%	abs(U) % todo

	if im
		im plc 1 2
		im(1, full(U))
		im(2, U.codes)
		drawnow
	end

	fatrix2_tests(U, 'complex', 1) % check complex data
	test_adjoint(U, 'complex', 1, 'tolre', 1e-9);
end

if 1 % 2d
	mask = true(8*3,16);
	mask(1:3) = false;
	level = 3;
	wname = 'haar';
	wname = 'sym2';
	U = Godwt1(mask, 'level', level, 'wname', wname);

	x = ellipse_im(size(mask));
	if im
		im plc 1 2
		im(1, x)
		im(2, U * x)
		drawnow
	end

	fatrix2_tests(U, 'complex', 1) % check complex data
	test_adjoint(U, 'complex', 1, 'tolre', 1e-9, 'big', 1);
end
