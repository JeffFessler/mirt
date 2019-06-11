 function ob = Gdown(Nd, varargin)
%function ob = Gdown(Nd, varargin)
%|
%| Construct Gdown object, which performs downsampling.
%| This is useful perhaps for iterative zooming methods.
%| See Gdown_test.m for example usage.
%|
%| in
%|	Nd	[1 D]		input signal dimensions
%|
%| options
%|	mask	logical			default: true(Nd)
%|	down	1|2|3|...		down sampling factor (default 2)
%|	type	func|Gsparse|...	not done.  default: 'func'
%| out
%|	ob			fatrix2 object
%|
%| Copyright 2006-8-25, Jeff Fessler, University of Michigan

if nargin == 1 && streq(Nd, 'test'), Gdown_test, return, end
if nargin < 1, ir_usage, end

% defaults
arg.mask = [];
arg.idim = Nd;
arg.down = 2;
arg.type = 'func';
arg.class = 'fatrix2';

% options
arg = vararg_pair(arg, varargin);

arg.odim = arg.idim ./ arg.down;
tmp = round(arg.odim) * arg.down;
if any(tmp ~= arg.idim)
	error 'only integer multiple image size supported'
end

if length(arg.idim) > 2
	error 'only 2d implemented due to downsample2() limitation'
end

if isempty(arg.mask)
	arg.mask = true(arg.idim);
end

arg.up_scale = 1 / arg.down^numel(arg.idim);

switch arg.class
case 'Fatrix'
	arg.np = sum(arg.mask(:));
	arg.dim = [prod(arg.odim) arg.np]; % nd x np
	ob = Fatrix(arg.dim, arg, ...
		'abs', @(ob) ob, ...
		'forw', @Gdown_forw_Fatrix, 'back', @Gdown_back_Fatrix);
case 'fatrix2'
	forw = @(arg, x) downsample2(x, arg.down);
	back = @(arg, y) fatrix2_maskit(arg.mask, ...
			upsample_rep(arg.up_scale * y, arg.down));
	ob = fatrix2('arg', arg, 'forw', forw, 'back', back, ...
		'abs', @(ob) ob);
		
otherwise
	fail 'not done'
end



% Gdown_forw_Fatrix(): y = A * x
% in
%	x	[np L] or [(Nd) L]
% out
%	y	[M L]
%
function y = Gdown_forw_Fatrix(arg, x)

[x ei] = embed_in(x, arg.mask, arg.np);

LL = size(x, 1+length(arg.idim)); % *L
if LL == 1
	y = downsample2(x, arg.down);
else
	y = [];
	for ll=1:LL
		y(:,:,ll) = downsample2(x(:,:,ll), arg.down); % fix: generalize
	end
end

y = ei.shape(y);


% Gdown_back_Fatrix(): x = A' * y
% in
%	y	[M L]
% out
%	x	[np L]
%
function x = Gdown_back_Fatrix(arg, y)

[y eo] = embed_out(y, arg.odim);

LL = size(y, 1+length(arg.odim)); % *L
if LL == 1
	x = upsample_rep(y, arg.down);
else
	for ll=1:LL
		x(:,:,ll) = upsample_rep(y(:,:,ll), arg.down); % fix: generalize
	end
end
x = x * arg.up_scale;

x = eo.shape(x, arg.mask, arg.np);


%
% Gdown_test
%
function Gdown_test

dim = [4 6];
down = 2;
A = Gdown(dim, 'down', down);

x = reshape(1:prod(dim), dim);
y1 = downsample2(x, down);
y2 = A * x;
jf_equal(y1, y2)

y = y1;
x1 = upsample_rep(y, down) * A.arg.up_scale;
x2 = A' * y;
jf_equal(x1, x2)

mask = true(dim);
if isa(A, 'Fatrix')
	Fatrix_test_basic(A, mask)
else
	fatrix2_tests(A)
end

test_adjoint(A);
