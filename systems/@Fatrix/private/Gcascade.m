  function ob = Gcascade(arg1, arg2, varargin)
%|function ob = Gcascade(A1, A2, options)
%|
%| Construct Gcascade object, which is the cascade of two objects: A = A1 * A2.
%|
%| See Gcascade_test.m for example usage.
%|
%| in
%|	arg1	matrix | Fatrix		object that can do "mtimes" and "size"
%|	arg2	matrix | Fatrix		''
%|
%| options
%|	'chat'				verbosity
%|
%|	'blockify_data'			see Fatrix.m about block handlers
%|	'block_setup'			provide routine to handle Gblock(A)
%|	'mtimes_block'			provide routine to handle A{m}
%|					(built-in already for scalar * Fatrix)
%|
%| out
%|	ob	A1 * A2
%|
%| Copyright 2005-7-21, Jeff Fessler, University of Michigan

if nargin == 1 && streq(arg1, 'test'), Gcascade_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% defaults
arg.chat = 0;
arg.blockify_data = [];
arg.block_setup = [];
arg.mtimes_block = [];
arg = vararg_pair(arg, varargin);

arg.A1 = arg1;
arg.A2 = arg2;

if isnumeric(arg1) && isscalar(arg1) % scalar * object
	dim = size(arg.A2);
	if isfield(struct(arg.A2), 'arg') ...
	&& isfield(struct(arg.A2).arg, 'odim')
		arg.odim = arg.A2.arg.odim; % to enable Gblock
	end

	if isfield(struct(arg.A2), 'handle_mtimes_block')
		if isempty(arg.mtimes_block)
			arg.mtimes_block = @Gblock_scalar_mtimes_block;
		else
			fail 'option over-ride of mtimes_block?'
		end
	end

	if isfield(struct(arg.A2), 'handle_blockify_data')
		if isempty(arg.blockify_data)
			arg.blockify_data = arg.A2.handle_blockify_data;
		else
			fail 'option over-ride of blockify_data?'
		end
	end

	% todo: block_setup

else % A1 * A2
	if size(arg.A1, 2) ~= size(arg.A2, 1)
		error 'size mismatch'
	end
	dim = [size(arg.A1,1) size(arg.A2,2)];
end

%
% build Fatrix object
%
ob = Fatrix(dim, arg, 'caller', mfilename, ...
	'forw', @Gcascade_forw, 'back', @Gcascade_back, ...
	'power', @Gcascade_power, ...
	'block_setup', arg.block_setup, ...
	'mtimes_block', arg.mtimes_block);


%
% Gcascade_forw(): y = A * x
%
function y = Gcascade_forw(arg, x)

y = arg.A1 * (arg.A2 * x);


%
% Gcascade_back(): x = A' * y
%
function x = Gcascade_back(arg, y)

x = arg.A2' * (arg.A1' * y);


%
% Gcascade_power(): A.^p
%
function ob = Gcascade_power(ob, pow)
if isnumeric(ob.arg.A1) && isscalar(ob.arg.A1)
	ob.arg.A1 = ob.arg.A1 .^ pow;
	ob.arg.A2 = ob.arg.A2 .^ pow;
else
	error 'power defined only for cascade of scalar * object'
end


%
% Gblock_scalar_mtimes_block()
%
function out = Gblock_scalar_mtimes_block(arg, is_transpose, x, iblock, nblock)
out = arg.A2.handle_mtimes_block(arg.A2.arg, is_transpose, x, iblock, nblock);
out = arg.A1 * out; % multiply by the scalar


%
% Gcascade_test
%
function Gcascade_test
a = [1:4]';
b = [2:5]';
c = diag_sp(a);
d = diag_sp(b);
e = Gcascade(c, d);
%e .^ 2;
f = [3:6]';
equivs(e * f, a .* b .* f)
g = Gcascade(7, c); % scalar * object
equivs(g * f, 7 * a .* f)

h = c * d; % check Fatrix * Fatrix
equivs(h * f, e * f)

equivs(g.^3 * f, 7^3 * a.^3 .* f) % power
