 function ob = Gmatrix(matrix, varargin)
%function ob = Gmatrix(matrix, [options])
%|
%| Construct fatrix2 object from an object (such as a matrix) that can do
%| multiplication and transpose.   The point of this would be to endow
%| the matrix with the masking features provided by fatrix2.
%|
%| in
%|	matrix	[M N]	could be a matrix or a Fatrix for example.
%|
%| option		(see fatrix2 for descriptions)
%|	idim	[]	default: N
%|	imask	[]	default true(N,1)
%|	odim	[]	default: M
%|	omask	[]	default: true(M,1)
%|
%| out
%|	ob	fatrix2: [M N]
%|
%| Copyright 2010-12-23, Jeff Fessler, University of Michigan

if nargin == 1 && streq(matrix, 'test'), Gmatrix_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.idim = [];
arg.imask = [];
arg.odim = [];
arg.omask = [];

arg = vararg_pair(arg, varargin);

[M N] = size(matrix);

if isempty(arg.idim) && isempty(arg.imask)
	arg.idim = N;
end

if isempty(arg.odim) && isempty(arg.omask)
	arg.odim = M;
end

if isempty(arg.imask)
	arg.apply_imask = @(arg, x) col(x);
else
	arg.apply_imask = @(arg, x) x(arg.imask);
end

if isempty(arg.imask)
	if N ~= prod(arg.idim)
		fail('idim mismatch: N=%d idim=[%s]', N, mat2str(arg.idim))
	end
else
	if N ~= sum(arg.imask(:))
		fail('idim mismatch: N=%d sum(imask)=%d', N, sum(arg.imask(:)))
	end
end

if isempty(arg.omask)
	if M ~= prod(arg.odim), fail('odim mismatch'), end
else
	if M ~= sum(arg.omask(:)), fail('odim mismatch'), end
end

arg.matrix = matrix;

if issparse(matrix) % stupid matlab requires sparse * double
	forw = @(arg,x) fatrix2_embed(arg.omask, arg.odim, ...
		arg.matrix * double(fatrix2_masker(arg.imask, arg.idim, x)));
	back = @(arg,y) fatrix2_embed(arg.imask, arg.idim, ...
		arg.matrix' * double(fatrix2_masker(arg.omask, arg.odim, y)));
else
	forw = @(arg,x) fatrix2_embed(arg.omask, arg.odim, arg.matrix ...
		* fatrix2_masker(arg.imask, arg.idim, x));
	back = @(arg,y) fatrix2_embed(arg.imask, arg.idim, arg.matrix' ...
		* fatrix2_masker(arg.omask, arg.odim, y));
end

ob = fatrix2('arg', arg, ...
	'idim', arg.idim, 'imask', arg.imask, ...
	'odim', arg.odim, 'omask', arg.omask, ...
	'forw', forw, 'back', back, ...
	'forw_block', @Gmatrix_forw_block, ...
	'back_block',  @Gmatrix_back_block, ...
	'abs', @Gmatrix_abs, 'power', @Gmatrix_power);


% generalization of masker for empty mask
% [(dim)] -> [np]
function z = fatrix2_masker(mask, dim, x)
if isempty(mask)
	z = x(:);
else
	z = x(mask);
end


% Gmatrix_abs()
% |A|
function ob = Gmatrix_abs(ob, p)
ob.arg.matrix = abs(ob.arg.matrix);


% Gmatrix_power()
% A .^ p
function ob = Gmatrix_power(ob, p)
ob.arg.matrix = ob.arg.matrix .^ p;


% Gmatrix_forw_block()
function y = Gmatrix_forw_block(arg, x, iblock, nblock)

x = fatrix2_masker(arg.imask, arg.idim, x);
if issparse(arg.matrix) % stupid matlab requires sparse * double
	x = double(x);
end

na = arg.odim(end); % last dimension
n1 = prod(arg.odim(1:end-1)); % first dimension(s)
ia = iblock:nblock:na;
[ii1 iia] = ndgrid(1:n1, ia);
ii = sub2ind([n1 na], ii1, iia);
tmp = arg.matrix(ii,:);
y = tmp * x(:);

if ~isempty(arg.omask)
	fail 'both omask and block not supported'
end
odim = arg.odim; odim(end) = numel(ia);
y = fatrix2_embed(arg.omask, odim, y);


% Gmatrix_back_block()
% expects y to be 'compact' size with size(y,end) = na
function x = Gmatrix_back_block(arg, y, iblock, nblock)

odim = []; % trick: fatrix2_masker() ignores odim
y = fatrix2_masker(arg.omask, odim, y);
if issparse(arg.matrix) % stupid matlab requires sparse * double
	y = double(y);
end

na = arg.odim(end); % last dimension
n1 = prod(arg.odim(1:end-1)); % first dimension(s)
ia = iblock:nblock:na;
[i1 ia] = ndgrid(1:n1, ia);
ii = sub2ind([n1 na], i1, ia);
tmp = arg.matrix(ii,:);
x = tmp' * y(:);
x = fatrix2_embed(arg.imask, arg.idim, x);


%
% Gmatrix_test()
%
function Gmatrix_test

M = 5;
N = 4;
x = [1:N]' * 10;
x = single(x); % stress test sparse case

for ii=1:2
	if ii == 1
		A = reshape(1:N*M, [M N]);
	else
		A = sparse(A);
	end

	if 1
		B = Gmatrix(A); % no mask
		y1 = A * double(x);
		y2 = B * x;
		jf_equal(y1, y2)

		fatrix2_tests(B)
		test_adjoint(B);
	end

	if 1
		imask = true(3,2); imask(1:2) = false;
		omask = true(2,3); omask(1) = false;
		B = Gmatrix(A, 'imask', imask, 'omask', omask); % masked

		y1 = A * double(x);
		y2 = B * x;
		jf_equal(y1, y2)

		fatrix2_tests(B, 'halt', 0)
		test_adjoint(B);
	end
end
