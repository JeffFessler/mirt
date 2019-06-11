 function ob = Gtranslate(mask, varargin)
%function ob = Gtranslate(mask, options)
%|
%| Construct Gtranslate object for image registration.
%| (Internally stores image-sized arrays so not very memory efficient.)
%|
%| See Gtranslate_test() below for example usage.
%|
%| in
%|	mask	size(image)	logical array of object support.
%|
%| options
%|	'shift'	[ndim]		shift amount (default 0)
%|	'type'	char		translation method:
%|				'circshift' (integer shifts only) (default)
%|				'fft' (terrible for non-integer shifts)
%|				'interpn,linear,circ' linear interpolation
%|				with circulant boundary conditions
%|
%|	these methods do not have a properly matched adjoint (todo):
%|				'bspline3'
%|				'interpn,linear' linear interpolation
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|			nd = np for 'conv,same' type
%|
%| Copyright 2010-3-22, Jeff Fessler, University of Michigan
%|
%| 2012-8-1 'interpn,linear,circ' option added by F. Zhao and M. Muckley

if nargin == 1 && streq(mask, 'test'), Gtranslate_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = mask;

% option defaults
arg.shift = 0; % identity matrix
arg.type = 'circshift';
arg.class = 'fatrix2';

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

arg.shift = arg.shift(:)'; % row vector

arg.ndim = ndims(mask);
if arg.ndim == 2 && size(mask,2) == 1
	arg.ndim = 1;
end
if numel(arg.shift) ~= arg.ndim
	fail 'ndim mismatch'
end

switch arg.type

case 'bspline3'
	if ndims(mask) ~= 2, fail 'only 2d done', end

	ig = image_geom('nx', size(mask,1), 'dx', 1, ...
			'ny', size(mask,2), 'dy', 1);
	kg = knot_geom( 'nx', 1, 'mx', 8*ig.nx, 'offset_x', ig.nx/2, ...
			'ny', 1, 'my', 8*ig.ny, 'offset_y', ig.ny/2);
	Bx = makeB(ig, kg);
	By = Bx;
	alphax = zeros(kg.dim, 'single');
	alphay = zeros(kg.dim, 'single');
	alphax(1) = -2.25 * arg.shift(1); % trick: empirical
	if arg.ndim > 1
		alphay(1) = -2.25 * arg.shift(2);
	else
		alphay(1) = 0; % 1d
	end
	arg.W = makeW({Bx, By}, {alphax, alphay});

	arg.fun_forw = @(arg, x) ...
		reshape(arg.W * BsplVal2CoMirr(single(x)), size(x));
	arg.fun_back = @(arg, y) arg.mask .* ...
		reshape(arg.W' * y, size(y)); % todo: wrong
%	arg.fun_back = @(arg, y) BsplCo2ValTranMirr(arg.W' * y);
%	arg.fun_back = @(arg, y) fail 'bspline3 adjoint not done';
	warn 'bspline3 adjoint not done to match!';

case 'circshift'
	if any(round(arg.shift) ~= arg.shift) % need to use fft?
		fail 'circshift needs integer shifts'
	end
	arg.fun_forw = @(arg, x) circshift(x, arg.shift);
	arg.fun_back = @(arg, y) arg.mask .* circshift(y, -arg.shift);

case 'fft'
	Nd = size(mask);
	for id = 1:arg.ndim
		N = Nd(id);
		kk = single([0:N-1]);
		phase{id} = exp(-2i * pi / N * kk * arg.shift(id));
	end
	tmp = ndgrid_jf('mat', phase);
	tmp = prod(tmp, 1 + arg.ndim);
	arg.phase = single(tmp);

%	arg.fun_forw = @(x) ifftn(fftn(x) .* phase);
	arg.fun_forw = @(arg, x) translate_fft(x, arg.phase);
	arg.fun_back = @(arg, y) arg.mask .* translate_fft(y, conj(arg.phase));

case 'interpn,linear'
	Nd = size(mask);
	for id = 1:arg.ndim
		N = Nd(id);
		xf{id} = single([1:N]) - arg.shift(id);
		xb{id} = single([1:N]) + arg.shift(id);
	end
	if ndims(mask) == 2 && size(mask,2) == 1 % 1d
		arg.xf = col(xf{1});
		arg.xb = col(xb{1});
		arg.fun_forw = @(arg, x) interp1(x, arg.xf, 'linear', 0);
		arg.fun_back = @(arg, y) interp1(y, arg.xb, 'linear', 0);
	else
		arg.xf = ndgrid_jf('cell', xf);
		arg.xb = ndgrid_jf('cell', xb);
		arg.fun_forw = @(arg, x) interpn(x, arg.xf{:}, 'linear', 0);
		arg.fun_back = @(arg, y) arg.mask .* ...
			interpn(y, arg.xb{:}, 'linear', 0); % todo: wrong
	end
	warn 'interpn adjoint not done to match!';
%	arg.fun_back = @(arg, y) fail 'interp adjoint not done';

case 'interpn,linear,circ'
	Nd = size(mask) + 2 * ceil(abs(arg.shift));
	for id = 1:arg.ndim
		N = Nd(id);
		xf{id} = single([1:N]) - arg.shift(id);
		xb{id} = single([1:N]) + arg.shift(id);
	end
	arg.xf = ndgrid_jf('cell', xf);
	arg.xb = ndgrid_jf('cell', xb);
	arg.fun_forw = @(arg, x) translate_interp_circ_forw(arg, x);
	arg.fun_back = @(arg, y) arg.mask .* translate_interp_circ_back(arg, y);

otherwise
	fail('type "%s" unknown', arg.type)
end


% build object
switch arg.class

case 'Fatrix'
	arg.np = sum(mask(:));
	arg.nd = prod(size(mask));
	dim = [arg.nd arg.np];
	ob = Fatrix(dim, arg, 'caller', 'Gtranslate', ...
		'forw', @Gtranslate_forw_Fatrix, ...
		'back', @Gtranslate_back_Fatrix);

case 'fatrix2'
	idim = size(mask);
	if numel(idim) == 2 && idim(2) == 1
		idim = idim(1); % 1d
	end
	ob = fatrix2('mask', mask, 'arg', arg, ...
		'idim', idim, 'odim', idim, ...
		'forw', arg.fun_forw, 'back', arg.fun_back);

otherwise
	fail('class "%s" unknown', arg.class)
end


% Gtranslate_forw_Fatrix()
% y = A * x
function y = Gtranslate_forw_Fatrix(arg, x)

[x ei] = embed_in(x, arg.mask, arg.np);

if ndims(x) > arg.ndim
	y = zeros(arg.nd, 1, class(x));
	for ll=1:size(x, ndims(x))
		y(:,ll) = col(arg.fun_forw(arg, stackpick(x, ll)));
	end
	if ndims(arg.mask) > 2 || size(arg.mask, ndims(arg.mask)) > 1 % > 1d
		y = reshapee(y, size(arg.mask), []);
	end
else
	y = arg.fun_forw(arg, x);
end

y = ei.shape(y);


% Gtranslate_back_Fatrix()
% x = A' * y (adjoint)
function x = Gtranslate_back_Fatrix(arg, y)

[y eo] = embed_out(y, size(arg.mask));

if ndims(y) > ndims(arg.mask)
	x = zeros(arg.nd, 1, class(y)); % trick!
	for ll=1:size(y, ndims(y))
		x(:,ll) = col(arg.fun_back(arg, stackpick(y, ll)));
	end
	x = reshapee(x, size(arg.mask), []);
else
	x = arg.fun_back(arg, y);
end

x = eo.shape(x, arg.mask, arg.np);


% translate_fft()
function y = translate_fft(x, phase)
y = ifftn(fftn(x) .* phase);
if isreal(x)
	y = real(y);
end


% translate_interp_circ_forw()
function y = translate_interp_circ_forw(arg, x)
padsize = ceil(abs(arg.shift));
%x1 = padarray(x, padsize, 'circular', 'both'); % needs image toolbox
if numel(arg.shift) == 1
	npad = 2*padsize + size(x,1);
else
	npad = 2*padsize + size(x);
end
x1 = ir_pad_into_center(x, npad, 'circ', 1);
if numel(arg.shift) == 1
	y = interp1(x1, arg.xf{1}, 'linear', 0);
else
	y = interpn(x1, arg.xf{:}, 'linear', 0);
end
%msk = logical(padarray(true(size(x)), padsize)); % image toolbox
msk = ir_pad_into_center(true(size(x)), npad);
y = reshape(y(msk), size(x));


% translate_interp_circ_back()
function x = translate_interp_circ_back(arg, y)
padsize = ceil(abs(arg.shift));
%y1 = padarray(y, padsize, 'circular', 'both'); % needs image toolbox
if numel(arg.shift) == 1
	npad = 2*padsize + size(y,1);
else
	npad = 2*padsize + size(y);
end
y1 = ir_pad_into_center(y, npad, 'circ', 1);
if numel(arg.shift) == 1
	x = interp1(y1, arg.xb{1}, 'linear', 0);
else
	x = interpn(y1, arg.xb{:}, 'linear', 0);
end
%msk = logical(padarray(true(size(y)), padsize)); % image toolbox
msk = ir_pad_into_center(true(size(y)), npad);
x = reshape(x(msk), size(y));


% Gtranslate_test()
function Gtranslate_test

rng(0)
classes = {'fatrix2', 'Fatrix'};
for ic=1:2
for id=1:2
	if id == 1
		mask = true(8,1); ishift = 2; fshift = [2.3]; % 1d
	else
		mask = true(10,8); ishift = [7 3]; fshift = [2.3 1.7]; % 2d
	end
	mask(1) = false;
	args = {mask, 'class', classes{ic}};

	x = mask .* rand(size(mask));

	btypes = {'bspline3', 'interpn,linear'}; % non-matched adjoint cases
	ftypes = {'fft', 'interpn,linear,circ'}; % arbitrary shifts with adjoint
	itypes = {'circshift'}; % integer shifts (orthonormal operation)

	for it = 1:numel(btypes)
		btype = btypes{it};
		A = Gtranslate(args{:}, 'type', btype, 'shift', fshift);
		if 0 % todo
			fatrix2_tests(A, 'complex', 0, 'halt', 0, ...
			'check1', false, 'full', false) % because of bad adjoint
		end

		if 0
			a1 = full(A)
			a2 = full(A')'
		end
%		test_adjoint(A, 'complex', 1); % todo: fails due to bad adjoint
	end

	for it = 1:numel(ftypes)
		ftype = ftypes{it};
		A = Gtranslate(args{:}, 'type', ftype, 'shift', fshift);
		fatrix2_tests(A, 'complex', 1)
		test_adjoint(A, 'complex', 1);
	end

	for it = 1:numel(itypes)
		itype = itypes{it};
		A = Gtranslate(args{:}, 'type', itype, 'shift', ishift);

		if id > 1
			jf_equal(x, A' * (A * x)) % orthonormal
		end
		jf_equal(x(mask), A' * (A * x(mask))) % orthonormal
		fatrix2_tests(A, 'complex', 1)
		test_adjoint(A, 'complex', 1);
	end
end
end

if 0
	tmp = 0 * mask; tmp(ceil(end/2), ceil(end/2)) = 1;
	im plc 2 3
	im(1, tmp)
	tmp = A * tmp;
	im(2, abs(tmp), 'abs'), cbar
	im(3, real(tmp), 'real'), cbar
	im(4, imag(tmp), 'imag'), cbar
return
end

% try it for image registration
if 1
%	f.type = 'fft'; % bad results!
%	f.type = 'bspline3';
%	f.type = 'interpn,linear';
	f.type = 'interpn,linear,circ';

	f1 = ellipse_im(128, [], 'oversample', 2);

	im plc 2 2
	im(1, f1), cbar

	mask = true(size(f1));
	A = Gtranslate(mask, 'type', f.type, 'shift', fshift);
	y1 = A * f1;
	im(2, y1), cbar

	if 0 % examine integer shifts
		Ac = Gtranslate(mask, 'type', 'circshift', 'shift', ishift);
		yc = Ac * f1;
		Ab = Gtranslate(mask, 'type', f.type, 'shift', ishift);
		yb = Ab * f1;
		yb = reshape(yb, size(mask));
		im clf, im_toggle(yb, yc, [0 255])
	return
	end

	if 0 % examine shifts
		Af = Gtranslate(mask, 'type', f.type, 'shift', fshift);
		yf = Af * f1;
		Ab = Gtranslate(mask, 'type', 'bspline3', 'shift', fshift);
		yb = Ab * f1;
		yb = reshape(yb, size(mask));
		im clf, im_toggle(yf, yb, [0 255])
	return
	end

	shifts = linspace(0, 4, 21);
	cost = zeros(size(shifts));
	for is=1:length(shifts)
		tmp = [shifts(is) fshift(2)]; % in dim1 only
		A2 = Gtranslate(mask, 'type', f.type, 'shift', tmp);
		y2 = A2 * f1;
		cost(is) = norm(y2(:) - y1(:))^2;
	end

	if im
		im subplot 3
		plot(shifts, cost, '-o', fshift(1), 0, 'rx')
	end
end
