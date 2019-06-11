  function ob = Gwave2(varargin)
%|function ob = Gwave2([args])
%|
%| Construct object that performs wavelet analysis of a 2D image,
%| for regularization based on sparsity such as in compressed sensing.
%|
%| option
%|	'mask'	logical [nx ny]	image-domain mask, often: true(nx,ny)
%|	'wname'		wavelet name (default: 'haar')
%|				(anything else requires wavelet toolbox) 
%|	'nlevel' 	integer specifying decomposition level (default 2)
%|	'scale_filters' 	factor for scaling wavelet filters
%|				(default: 1/sqrt(2))
%|	'redundancy' 	char	default: 'undecimated' (only option for now)
%|
%| out:
%|	ob	[np]	fatrix2 object, where np = sum(mask(:))
%|
%| Copyright 2011-01-19, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), Gwave2_test, return, end

% options
arg.mask = [];
arg.dwtmode = 'per'; % periodic boundaries
arg.redundancy = 'undecimated'; % undecimated wavelet transform using wavelet filters corresponding to standard orthonormal wavelets
arg.nlevel = 2; % wavelet transform with 2 levels 
arg.includeApprox = false; % do not include approximation level in the regularizer
arg.wname = 'haar';
arg.scale_filters = 1 / sqrt(2); % avoids a product with 0.5 during inverse undecimated wavelet transform 

arg = vararg_pair(arg, varargin);

if isempty(arg.mask), fail 'must provide a mask', end
if ~islogical(arg.mask), error 'mask must be logical', end

arg = Gwave2_setup1(arg);
%filt = Gwave2_setup(arg.wname, arg.scale_filters)


% build fatrix2 object
ob = fatrix2('imask', arg.mask, 'odim', arg.odim, 'arg', arg, ...
	'abs', @Gwave2_abs, 'power', @Gwave2_power, ...
	'forw', @Gwave2_forw, 'back', @Gwave2_back);


% Gwave2_setup1()
% version based on sathish ramani's code
function arg = Gwave2_setup1(arg)

dwtmode(arg.dwtmode, 'nodisp');

[filt.lod, filt.hid, filt.lor, filt.hir] = wfilters(arg.wname); % filters

filt.lod = arg.scale_filters * filt.lod;
filt.hid = arg.scale_filters * filt.hid;
filt.lor = arg.scale_filters * filt.lor;
filt.hir = arg.scale_filters * filt.hir;

arg.filt = filt;

switch arg.redundancy
case 'undecimated'
	[nx ny] = size(arg.mask);
	if arg.includeApprox
		arg.odim = [nx ny 3*arg.nlevel+1];
	else
		arg.odim = [nx ny 3*arg.nlevel];
	end
otherwise
	fail('redundancy %s not done', arg.redunancy)
end


%{
% Gwave2_setup()
% future version
function filt = Gwave2_setup(wname, scale_filters)

if exist('dwtmode') ~= 2
	fail 'Gwave2 requires matlab wavelet toolbox'
return
end

if streq(wname, 'haar')
	filt.lod = [1 1] / sqrt(2);
	filt.hid = [-1 1] / sqrt(2);
	filt.lor = [1 1] / sqrt(2);
	filt.hir = [1 -1] / sqrt(2);
else
	[filt.lod, filt.hid, filt.lor, filt.hir] = wfilters(wname);
end

filt.lod = scale_filters * filt.lod;
filt.hid = scale_filters * filt.hid;
filt.lor = scale_filters * filt.lor;
filt.hir = scale_filters * filt.hir;

%}



% Gwave2_forw(): y = A * x 
% wavelet analysis on x.
function y = Gwave2_forw(arg, x)

switch arg.redundancy
case 'undecimated'
	y = myswt2(x, arg.nlevel, arg.filt.lod, arg.filt.hid);
	if ~arg.includeApprox
		y = y(:,:,1:3*arg.nlevel);
	end
otherwise
	fail('redundancy %s not done', arg.redundancy)
	if ~isreal(x)
		fail 'todo'
	end
end


% Gwave_back(): x = A' * y
% adjoint of wavelet analysis
function x = Gwave2_back(arg, y)

switch arg.redundancy
case 'undecimated'
	if ~arg.includeApprox
		[nx ny] = size(arg.mask);
		y(:, :, end+1) = zeros(nx, ny);
	end
	x = myiswt2(y, arg.filt.lor, arg.filt.hir);
otherwise
	fail('redundancy %s not done', arg.redundancy)
	if ~isreal(x)
		fail 'todo'
	end
end
x = x .* arg.mask;


% Gwave2_abs()
% for C = abs(C)
function ob = Gwave2_abs(ob)
fun = @(x) abs(x);
filt.lod = fun(ob.arg.filt.lod);
filt.hid = fun(ob.arg.filt.hid);
filt.lor = fun(ob.arg.filt.lor);
filt.hir = fun(ob.arg.filt.hir);
ob.arg.filt = filt;


% Gwave2_power()
% for C = C.^p
function ob = Gwave2_power(ob, p)
fun = @(x) x .^ p;
filt.lod = fun(ob.arg.filt.lod);
filt.hid = fun(ob.arg.filt.hid);
filt.lor = fun(ob.arg.filt.lor);
filt.hir = fun(ob.arg.filt.hir);
ob.arg.filt = filt;


% Gwave2_test()
function Gwave2_test

if 1 % look at tiny case
	mask = true(2^3+2,2^2);
	C = Gwave2('mask', mask);
	Cf = full(C);
	tmp = embed(Cf', mask);
	im(tmp)
%	jf_slicer(tmp, 'iz', 1);
prompt
end

mask = true(2^5,2^4); mask(1) = false; % stress
C = Gwave2('mask', mask);

x = zeros(size(mask)); x(end/2+1,floor(end/2+1)) = 1;
tmp = C * x;
im(tmp), cbar

% todo: test complex case too

if 1
	fatrix2_tests(C)
	test_adjoint(C, 'big', 1)
end

Ca = abs(C);
fatrix2_tests(Ca)
test_adjoint(Ca, 'big', 1)

C2 = C .^ 2;
fatrix2_tests(C2)
test_adjoint(C2, 'big', 1)
