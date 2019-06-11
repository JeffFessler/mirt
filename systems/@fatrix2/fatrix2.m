 function ob = fatrix2(varargin)
%function ob = fatrix2('odim', odim, 'mask', mask, 'arg', arg, [options])
%|
%| Construct fatrix2 object, a matrix generalization for representing
%| any linear operator.  The 'f' might stand for 'fake' or 'function-based'
%| or 'fancy' or maybe a last name?
%| The caller can provide a variety of (overloaded) "methods" for this object,
%| particularly multiplication and transposed multiplication.
%| Using this "meta object" allows the object designer to focus on the key
%| methods, rather than reinventing matrix methods for each new object.
%| See G*.m .m for examples.
%|
%| Basic methods include the following:
%|
%|	A * x		matrix-vector multiplication (returns a vector)
%|	A' * x		transposed multiplication
%|	A(:,7)		returns a column vector
%|	A(:,[7 9])	returns a fatrix2 with smaller size
%|	A(:,:)		returns a full (conventional) matrix
%|	full(A)		same as A(:,:)
%|	size(A)		[nd np]
%|	disp(A)
%|	sum(A)		1' * A
%|	A.method(args)	user-defined methods called as method(A, args)
%|
%|	The following operations return a new fatrix2 (if no errors):
%|	scalar * A
%|	A'		(Hermitian) transpose of A
%|	-A		uminus, same as (-1) * A
%|	A1 * A2		if A1 and A2 have compatible sizes
%|	A1 - A2		minus. if A1,A2 have same size
%|	A1 + A2		plus. if A1,A2 have same size
%|	[A1 A2 ...]	horzcat
%|	[A1; A2; ...]	vertcat
%|	A(:,1:7)	returns a new fatrix2 with size: [size(A,1) 7]
%|	A(1:9,:)	returns a new fatrix2 with size: [9 size(A,2)]
%|	abs(A)		if user-provides 'abs' function handle
%|	power(A)	A.^2 if user-provides 'power' function handle
%|	A{1}		block object (if supported via forw_block handle)
%|
%|	array-vs-vector conversion (largely for internal use):
%|	iembed(A, x)	fatrix2_embed(A.imask, A.idim, x)
%|	oselect(A, y)	fatrix2_select(A.omask, A.odim, y)
%|	(note that imask/omask and idim/odim are swapped on ctranspose)
%|
%| The operation y = ob * x works differently depending on the
%| dimensions of x, and the way ob was created, as follows.
%|
%| Let np = sum(mask(:)) and N = size(mask)
%| and let nd be the vector-mode output size, usually prod(odim).
%| If size(x) is [np 1] then y is [nd 1].
%|	This is the usual matrix-vector way (vector-mode).
%| If size(x) is [np (L)] then y is [nd (L)].  Obvious vector-mode extension.
%| If size(x) is [N (L)] then y is [odim (L)]. The "array mode"
%| This "y" will be zero outside of omask, if "forw()" works properly.
%| Caution: if such an "x" has nonzero values outside of imask then the
%| result may be unpredictable.
%|
%| Likewise, the operation x = ob' * y has multiple versions.
%| If size(y) is [nd 1] then x is [np 1].  The usual vector-mode way.
%| If size(y) is [nd (L)] then x is [np (L)].  Obvious vector-mode extension.
%| If size(y) is [odim (L)] then x is [N (L)]. The "array mode"
%| This "x" will be zero outside of imask, if "back()" works properly.
%| Caution: if such an "y" has nonzero values outside of omask then the
%| result may be unpredictable.
%|
%| The "array mode" is supported only if mask is truly 2d or higher
%| or for a 1d mask that is not all true.
%|
%| dimension "options" (some combination of which is required)
%|	'idim' []	input array dimensions, often just [np]
%|			if empty, then this will look for one inside arg.idim
%|			or if not there, will use size(imask)
%|	'odim' []	output array dimensions, often just [nd]
%|			if empty, then this will look for one inside arg.odim
%|	'[i]mask' []	logical array for input dimensions, e.g. true(np,1)
%|			if empty, then this will look for one inside arg.mask
%|			and if not there, will default to true(idim)
%|	'omask' []	logical array for ouput dimensions, e.g. true(nd,1)
%|			if empty, then this will default to true(odim)
%|
%| typical "options"
%|	'arg'		arguments passed to handle functions
%|				(usually a cell or struct)
%|	'meth'		{'name1', @handle1, 'doc1'; 'name2', @handle2, ...}
%|			other methods to be invoked via ob.name()
%|
%| handle "options" akin to standard matrices
%|	'forw'		forw(arg, x)	ob * x
%|			input: size(x) == idim.  size(output) == odim.
%|			output must be zero outside of omask (for array mode)
%|	'back'		back(arg, y)	ob' * y
%|			input: size(y) == odim.  size(output) == idim.
%|			output must be zero outside of imask (for array mode)
%|	'free'		free(arg)
%|	'abs'		ob = abs(ob): absolute value operation
%|	'power'		ob = power(ob, sup)	ob.^sup
%|	'sparse'	ob = sparse(ob): return sparse matrix (default loops)
%|
%| handle "options" for special capabilities
%|	'gram'		build_gram(ob, W, reuse), build A' * W * A
%|	'block_setup'	ob = block_setup(ob, varargin)
%|	'blockify_data'	cell_data = blockify_data(ob, array_data, varargin)
%|	'forw_block'	y = forw_block(arg, x, iblock, nblock)
%|	'back_block'	x = back_block(arg, y, iblock, nblock)
%|
%| special options
%|	'caller', str	name of calling routine ("meta class")
%|			default: name determined automatically by caller_name.
%|	'does_many'	if forw() method can do [idim *L] instead of just the
%|			usual [idim], and likewise for back() method.  default 0
%|
%| out
%|	ob		fatrix2 object
%|			special methods:
%|				ob.imask_array [idim] logical array
%|				ob.omask_array [odim] logical array
%|
%| Copyright 2004-6-29, Jeff Fessler, University of Michigan


% create default object, as required by Mathworks

% properties that the user can specify
ob.caller = '';
ob.arg = {};
%
ob.idim = [];
ob.odim = [];
ob.imask = [];
ob.omask = [];
%
ob.meth = struct([]);
ob.docs = struct([]);
ob.scale = 1; % scalar scale factor (optional)
ob.idiag = []; % [idim] diagonal scaling matrix (on right side)
ob.odiag = []; % [odim] diagonal scaling matrix (on left side)
ob.does_many = false;

% internal properties
ob.accept1d = false;
ob.size = [];
ob.nblock = [];

% trick: default to some simple functions that will work if
% "arg" is anything that knows how to multiply, e.g., a matrix.
% Must be a pointer to function (handle, not anonymous) so we can test
% if it is the default later, e.g., for full().
ob.handle_back = @fatrix2_def_back;
ob.handle_forw = @fatrix2_def_forw;
ob.handle_power = @fatrix2_def_pow;

ob.handle_abs = [];
ob.handle_free = [];
ob.handle_gram = []; % defer to default within build_gram.m
ob.handle_sparse = [];
ob.handle_forw_block = [];
ob.handle_back_block = [];
ob.handle_block_setup = [];
ob.handle_blockify_data = [];

if nargin == 0 % required by Mathworks
	if nargout == 0, help(mfilename), end
	ob = class(ob, 'fatrix2');
return
end

if nargin == 1 && isstruct(varargin{1}) % trick to convert struct to fatrix2
	ob = class(varargin{1}, 'fatrix2');
return
end

if nargin < 2 % because ('arg', arg) is required at least!
	help(mfilename)
	error(mfilename)
end

ob = vararg_pair(ob, varargin, 'subs', { ...
	'mask', 'imask';
	'abs', 'handle_abs';
	'free', 'handle_free';
	'back', 'handle_back';
	'forw', 'handle_forw';
	'gram', 'handle_gram';
	'power', 'handle_power';
	'sparse', 'handle_sparse';
	'forw_block', 'handle_forw_block';
	'back_block', 'handle_back_block';
	'block_setup', 'handle_block_setup';
	'blockify_data', 'handle_blockify_data'});

ob = fatrix2_setup_dim(ob);

if isempty(ob.caller)
	ob.caller = caller_name;
end

[ob.meth ob.docs] = fatrix2_methods(ob.meth, ob);

ob = class(ob, mfilename);
% ob.cleanup = onCleanup( @() free(ob) ); % todo: to invoke 'free'


% fatrix2_methods()
% (based on strum)
function [method, docs] = fatrix2_methods(methods, ob)

method = struct;
docs = {};

%method.imask_array = @fatrix2_imask_array;
%method.omask_array = @fatrix2_omask_array;

for ii=1:size(methods,1)
	name = methods{ii,1};
	handle = methods{ii,2};
	if isfield(ob, name) % methods names should differ from object fields
		fail('name conflict: %s', name)
	end
	if isfield(method, name)
		warn('name reuse: "%s"', name)
	end
	method.(name) = handle;

	doc = methods{ii,3};
	docs{end+1} = doc;
end
