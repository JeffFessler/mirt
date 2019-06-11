  function Fatrix_test_basic(A, mask, varargin)
%|function Fatrix_test_basic(A, mask, [option])
%|
%| perform basic tests of Fatrix objects
%| in
%|	A	[nd *N]		Fatrix object
%|	mask	[(N)]		logical support
%| option
%|	'x'	[(N)]		optional input test image
%|	'multi'	0|1		if 1, test multiple rhs: A * [x x]. default: 1
%|	'name'	char		name of input variable (optional for display)
%|	'caller' char		name of caller (optional for display)
%|	'halt'	0|1		if error, 1 (default) halt, 0 keyboard
%|	'complex' 0|1		if 1, test for complex inputs. default: 0
%|
%| Copyright 2005-8-2, Jeff Fessler, University of Michigan

if nargin == 1 && streq(A, 'test'), Fatrix_test_basic_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.x = [];
arg.multi = true;
arg.name = '';
arg.caller = '';
arg.complex = false;
arg.halt = true; % halt if error
arg = vararg_pair(arg, varargin);

x = arg.x;
if arg.complex && ~isempty(x) && isreal(x)
	fail 'asked for complex test but provided real x'
end

if isempty(x)
	switch ndims(mask)
	case 2
		[nx ny] = size(mask);
		ig = image_geom('nx', nx, 'ny', ny, 'dx', 1);
		x = ellipse_im(ig, 'shepplogan-emis', 'type', 'slow');
		clear ig nx ny

	case 3
		[nx ny nz] = size(mask);
		ig = image_geom('nx', nx, 'ny', ny, 'nz', nz, 'dx', 1, 'dz', 1);
		x = ellipsoid_im(ig, []);
		clear ig nx ny nz

	otherwise
		error 'only 2d and 3d implemented'
	end

	if arg.complex
		x = x + 1i * flipdim(x, 1); % something complex for testing
	end
	x = x .* mask; % test with mask-limited object
end

% sum
sum(sum(A));

% trick:
% some tests are meaningful only for objects that can go forward
% into an array instead of a vector, i.e., for tomography objects
% that make sinograms, but not for Gmri.

% array vs col
y1 = A * x;
flag_col_out = (size(y1,2) == 1); % for objects that always have col out
ydim = size(y1);
y2 = reshape(A * x(mask), ydim);
if arg.halt
	jf_equal(y1, y2)
else
	try
		jf_equal(y1, y2)
	catch
		warn 'A * x vs reshape(A * x(mask)) failed'
		keyboard
	end
end

if ~isreal(y1) && ~arg.complex
	warn 'recommend testing complex objects with ''complex'' option'
end

% A * [x() x()]
if arg.multi
	y2 = A * [x(mask) 2*x(mask)];
	z2 = [y1(:) 2*y1(:)];
	if arg.halt
		equivs(y2, z2) % may not be *exactly* equal due to 2*
	else
		try
			equivs(y2, z2)
		catch
			warn 'A * [x() x()] failed'
			keyboard
		end
	end
end

% version of ndims that returns "1" for column vectors!
ndims1 = @(x) ndims(x) - (size(x,ndims(x)) == 1);
catnext = @(a,b) cat(1+ndims1(a),a,b);

if isvar('A.arg.odim') % caution: handle single view projectors [ns nt 1 nrep]
	catp = @(a,b) cat(1+length(A.arg.odim), a, b);
else
	catp = catnext;
end

% A * [x x]
if arg.multi
	y2 = A * catnext(x, 2*x);
	z2 = catp(y1, 2*y1);
	if arg.halt
		equivs(y2, z2) % may not be *exactly* equal due to 2*
	else
		try
			equivs(y2, z2)
		catch
			warn 'A * [x x] failed'
			keyboard
		end
	end
end

y0 = y1;
if arg.complex && isreal(y0)
	warn 'complex testing is incomplete; tell jeff!'
end

% A' * y array vs col
x1 = A' * y0;
tmp = A' * y0(:);
x2 = embed(tmp, mask);
if flag_col_out % size(x1,1) == sum(mask(:)) % e.g., for Gmri
	x1 = embed(x1, mask);
end
if arg.halt
	jf_equal(x1, x2)
else
	try
		jf_equal(x1, x2)
	catch
		keyboard
	end
end

% A' * [y() y()]
if arg.multi
	x2 = A' * [y0(:) 2*y0(:)];
	v2 = [x1(mask) 2*x1(mask)];
	equivs(x2, v2) % may not be *exactly* equal due to 2*
end

% A' * [y y]
if arg.multi
	tmp = catp(y0, 2*y0);
	x2 = A' * catp(y0, 2*y0);
	v2 = catnext(x1, 2*x1);
	if flag_col_out
		x2 = embed(x2, mask);
	end
	equivs(x2, v2) % may not be *exactly* equal due to 2*
end

% indexing single and multiple arguments and test consistency thereof
% A(:,?)
c1 = A(:,1);
c2 = A(:,2);
c12 = A(:, [1 2]);
c12 = [c12(:,1) c12(:,2)]; % for fatrix2
jf_equal(c12, [c1 c2])

% A(?,:)
r1 = A(1,:);
r2 = A(2,:);
r12 = A([1 2], :);
r12 = [r12(1,:); r12(2,:)]; % for fatrix2
jf_equal(r12, [r1; r2])

if isempty(arg.name)
	arg.name = inputname(1);
	arg.caller = caller_name;
end
printm('passed for "%s" in %s', arg.name, arg.caller)


function Fatrix_test_basic_test

psf = [0 1 2 1 0; 1 2 4 3 1; 0 2 3 1 0];
psf = psf / sum(psf(:));
idim = [24 30];
mask = true(idim);
mask(1) = 0;
A = Gblur(mask, 'psf', psf); % todo: replace with a Fatrix!

Fatrix_test_basic(A, mask)
Af = A(:,:);
tmp = A(2,:);
jf_equal(tmp, Af(2,:))

A = Gblur(mask, 'psf', psf + 1i);
Fatrix_test_basic(A, mask, 'complex', 1)
