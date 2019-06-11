 function ob = Gdiag(diag, varargin)
%function ob = Gdiag(diag, [options])
%|
%| Construct diagonal "sparse" object, that does D * x = diag(d) * x via d .* x.
%|
%| in
%|	diag	[]	numeric array
%|
%| option
%|	'mask'	logical	with sum(mask(:)) = numel(diag)
%|	'class'		'fatrix2' (default) or 'Fatrix' (obsolete)
%|
%| out
%|	ob	fatrix2: [numel(diag) numel(diag)]
%|
%| Copyright 2005-4-6, Jeff Fessler, University of Michigan

if nargin == 1 && streq(diag, 'test'), Gdiag_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = [];
arg.class = 'fatrix2';

arg = vararg_pair(arg, varargin);

arg.diag = diag;
if ~isempty(arg.mask)
	arg.diag = embed(arg.diag, arg.mask); % [(N)]
end

arg.idim = size(arg.diag);
if arg.idim(end) == 1, arg.idim = arg.idim(1:end-1); end % 1d trick
arg.odim = arg.idim;

switch arg.class
case 'Fatrix'
	if isempty(arg.mask)
		arg.mask = true(size(arg.diag));
	end

	arg.odim = size(arg.mask); % Fatrix does not support omask
	arg.np = sum(arg.mask(:));
	arg.dim = [prod(arg.odim) arg.np];
	if any(~arg.mask(:))
		arg_gram = {}; % must use default internal gram
	else
		arg_gram = {'gram', @Gdiag_gram}; % can use simple version
	end
	ob = Fatrix(arg.dim, arg, ...
		'forw', @Gdiag_forw_Fatrix, 'back', @Gdiag_back_Fatrix, ...
		'abs', @Gdiag_abs, 'power', @Gdiag_power, arg_gram{:});

case 'fatrix2'

	forw = @(arg,x) bsxfun(@times, arg.diag, x);
	back = @(arg,y) bsxfun(@times, conj(arg.diag), y);

	block_arg = {};
	if numel(arg.idim) > 1 % only meaningful if >1D
		block_arg = {	'forw_block', @Gdiag_forw_block, ...
				'back_block', @Gdiag_back_block};
	end

	ob = fatrix2('arg', arg, 'does_many', 1, ...
		'imask', arg.mask, 'omask', arg.mask, ...
		'forw', forw, 'back', back, block_arg{:}, ...
		'abs', @Gdiag_abs, 'power', @Gdiag_power, ...
		'gram', @Gdiag_gram, 'sparse', @Gdiag_sparse);

otherwise
	fail 'class'
end


% Gdiag_forw_block()
function y = Gdiag_forw_block(arg, x, istart, nblock)
ia = istart:nblock:arg.odim(end);
d = reshape(arg.diag, [], arg.odim(end));
x = reshape(x, [], arg.odim(end));
d = d(:,ia);
x = x(:,ia);
tmp = d .* x;
y = reshape(tmp, [arg.odim(1:end-1) numel(ia)]);
% y = bsxfun(@times, y, x);


% Gdiag_back_block()
function x = Gdiag_back_block(arg, y, istart, nblock)
ia = istart:nblock:arg.odim(end);
d = reshape(arg.diag, [], arg.odim(end));
d = d(:,ia);
y = reshape(y, [], numel(ia));
tmp = conj(d) .* y;
x = zeros([prod(arg.idim(1:end-1)) arg.idim(end)], class(y));
x(:,ia) = tmp;
x = reshape(x, arg.idim);
%tmp = reshape(arg.diag, [], arg.odim(end));
%tmp = tmp(:,ia);
%tmp = reshape(tmp, [arg.odim(1:end-1) numel(ia)]);
%x = bsxfun(@times, conj(tmp), y);


% Gdiag_forw()
% for Fatrix
function y = Gdiag_forw_Fatrix(arg, x)

[x, ei] = embed_in(x, arg.mask, arg.np); % [(N) *L]
y = bsxfun(@times, arg.diag, x);
y = ei.shape(y);


% Gdiag_back_Fatrix(): x = A' * y
function x = Gdiag_back_Fatrix(arg, y)

[y, eo] = embed_out(y, arg.odim);
x = bsxfun(@times, conj(arg.diag), y);
x = eo.shape(x, arg.mask, arg.np);


% Gdiag_abs(): |D|
function ob = Gdiag_abs(ob)
ob.arg.diag = abs(ob.arg.diag);


% Gdiag_power(): D .^ p
function ob = Gdiag_power(ob, p)
ob.arg.diag = ob.arg.diag .^ p;


% Gdiag_sparse()
function sp = Gdiag_sparse(ob)
d = double(ob.arg.diag(:));
if ~isempty(ob.arg.mask)
	d = d(ob.arg.mask);
end
np = size(ob,2);
sp = sparse(1:np, 1:np, d, np, np, np);


% Gdiag_gram(): D'*W*D
% works only for fatrix2 or for Fatrix in case of all true mask
function [ob, reuse] = Gdiag_gram(ob, W, reuse, varargin)
if isempty(W)
	ob.arg.diag = abs(ob.arg.diag).^2;
else
	fail('W not implemented yet')
end


% Gdiag_test_block
function Gdiag_test_block
nb = 7;
na = 12;
d = reshape(1:(nb*na), nb, na) + 1i;
A = Gdiag(d, 'mask', []);
nblock = 3;
B = Gblock(A, nblock);
%pr size(B)
%pr size(B{2})
iblock = 2;
ia = iblock:nblock:na;
tmp1 = B{iblock} * col(ones(nb, na));
tmp2 = B{iblock} * ones(nb, na);
jf_equal(tmp1(:), tmp2(:))
jf_equal(tmp2, d(:,ia))

y = flipdim(d, 1);
tmp1 = B{iblock}' * col(y(:,ia));
tmp2 = B{iblock}' * y(:,ia);
jf_equal(tmp1(:), tmp2(:));
tmp3 = conj(d) .* y;
jf_equal(tmp2(:,ia), tmp3(:,ia))


% Gdiag_test()
function Gdiag_test

rng(0)
d = [-2:5]';
d = d * [1:3]; % 2d stress test
d = d + 1i*rand(size(d)); % complex
mask = true(size(d));
if 1 % mask stress
	mask(end-1) = false; % stress
	d = d(mask);
end

if 0 % test use of Gdiag as an "embed" object
	A = Gdiag(mask(mask), 'mask', mask);
	tmp1 = embed(d, mask); % [(N)]
	tmp2 = A * d; % [np]
	jf_equal(tmp1, tmp2) % fails
end

cl_list = {'Fatrix', 'fatrix2'};
for ic = 1:numel(cl_list)
	cl = cl_list{ic};

	D = Gdiag(d, 'mask', mask, 'class', cl);

	switch cl
	case 'Fatrix'
		Fatrix_test_basic(D, mask, 'complex', 1, 'halt', 0)
	case 'fatrix2'
		fatrix2_tests(D, 'complex', 1, 'halt', 0)
	end
	test_adjoint(D, 'complex', 1);

	x = d.^3;
if 0
	y = D*[x x]; y = y(:,1);
	jf_equal(d.*x, y)
	x = D'*[y y]; x = x(:,1);
	jf_equal(conj(d).*y, x)
	y = D.^2*x;
	jf_equal(d.^2.*x, y)
	y = abs(D)*x;
	jf_equal(abs(d).*x, y)
end

	T = build_gram(D, []);
	equivs(T * x, D' * (D * x));

	if 0 % scalar
		d = 7;
		D = Gdiag(d);
		x = ones(10,1);
		jf_equal(d * x, D * x);
	end
end

D1 = Gdiag(d+1, 'mask', mask);
D2 = Gdiag(d+2, 'mask', mask);
Dv = [D1; D2];
tmp1 = spdiag(d+1, 'nowarn');
tmp2 = spdiag(d+2, 'nowarn');
tmpv = sparse(Dv);
tmp0 = [tmp1; tmp2];
jf_equal(tmp0, tmpv)

Gdiag_test_block
