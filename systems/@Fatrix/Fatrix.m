  function ob = Fatrix(dim, arg, varargin)
%|function ob = Fatrix(dim, arg, handles, options)
%|
%| WARNING: any new objects should use fatrix2 instead of this!
%|
%| Construct Fatrix object, which is a 'fake matrix' or 'function-based matrix'
%| object, designed to generalize matrices by representing any linear operator.
%| The caller can provide a variety of (overloaded) "methods" for this object,
%| particularly multiplication and transposed multiplication.
%| Using this "meta object" allows the object designer to focus on the key
%| methods, rather than reinventing matrix methods for each new object.
%| Basic methods include things like "size" and "disp" etc., whereas more
%| advanced methods included ob(:,1) or [ob1; ob2] among many others.
%| See Gblur.m for examples.
%|
%|	ob * x		multiplication
%|	ob' * x		transposed multiplication
%|
%| It is up to the caller to provide function handles for both!
%|
%| in
%|	dim [2]		Fatrix "dimensions"
%|	arg		arguments passed to handle functions
%|				(usually a cell or struct)
%|
%| handles (all optional): 'name1', handle1, 'name2', handle2, ...
%|	'forw'		forw(arg, x)	ob * x
%|	'back'		back(arg, x)	ob' * x
%|	'gram'		build_gram(ob, W, reuse), build G' * W * G
%|	'free'		free(arg)
%|	'abs'		ob = abs(ob): absolute value operation
%|	'ufun'		out = ufun(ob, varargin) : user-defined function
%|	'block_setup'	ob = block_setup(ob, varargin)
%|	'blockify_data'	cell_data = blockify_data(ob, array_data, varargin)
%|	'mtimes_block'	mtimes_block(arg, is_transpose, x, iblock, nblock)
%|	'power'		ob = power(ob, sup)	ob.^sup
%|
%| options
%|	'caller', string	name of calling routine ("meta class")
%|				default: name determined automatically.
%|	'cascade_after', thing		ob * x -> thing * forw(arg, x)
%|	'cascade_before', thing		ob * x -> forw(arg, thing * x)
%|		these "things" can also be function handles (see mtimes_block.m)
%|
%| out
%|	ob		Fatrix object
%|
%| Copyright 2004-6-29, Jeff Fessler, University of Michigan


% create default object, as required by Mathworks
ob.caller = '';
ob.arg = {};
ob.dim = [];
ob.is_transpose = false;
ob.is_subref = false;
ob.nblock = [];
ob.iblock = [];

% trick: default to some simple anonymous functions that will work if
% "arg" is anything that knows how to multiply, e.g., a matrix.
ob.handle_back = @(M,y) M' * y;
ob.handle_forw = @(M,x) M * x;
ob.handle_power = @(M,p) M .^ p;

ob.handle_ufun = [];
ob.handle_abs = [];
ob.handle_free = [];
ob.handle_gram = [];
ob.handle_mtimes_block = [];
ob.handle_block_setup = [];
ob.handle_blockify_data = [];
ob.cascade_after = [];
ob.cascade_before = [];

if nargin == 0 % required by Mathworks
	if nargout == 0, help(mfilename), end
	ob = class(ob, 'Fatrix');
return
end

if nargin < 2
	help(mfilename)
	error(mfilename)
end

ob.arg = arg;
ob.dim = dim;

ob = vararg_pair(ob, varargin, 'subs', { ...
	'ufun', 'handle_ufun';
	'abs', 'handle_abs';
	'free', 'handle_free';
	'back', 'handle_back';
	'forw', 'handle_forw';
	'gram', 'handle_gram';
	'power', 'handle_power';
	'block_setup', 'handle_block_setup';
	'blockify_data', 'handle_blockify_data';
	'mtimes_block', 'handle_mtimes_block'});

if isempty(ob.caller)
	ob.caller = caller_name;
end

ob = class(ob, mfilename);
