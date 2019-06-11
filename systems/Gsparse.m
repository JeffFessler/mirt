 function ob = Gsparse(arg1, varargin)
%function ob = Gsparse(file.wtf | sparse | cell, options)
%|
%| Construct Gsparse object, either from a sparse matrix itself,
%| or from the arguments that would be passed to matlab's sparse() command,
%| or from an Aspire binary .wtf file.
%|
%| The Gsparse object overcomes some annoying limitations of
%| matlab's sparse() function and sparse datatype.  In particular, Gsparse
%| allows single precision arguments (that are converted silently to doubles,
%| as long as matlab continues to insist on that) instead of the useless
%| error message provided by Mathworks.  And you can multiply Gsparse object
%| times a non-double-precision vector (which is silently upgraded to doubles),
%| which Matlab does not support.
%|
%| More substantively, this object also supports the "multidimensional"
%| constructs needed in imaging problems. (support mask, subsets, etc.)
%|
%| This just uses an ordinary Matlab sparse matrix for the core!
%| See Gsparse_test.m for example usage.
%|
%| You create an system object by calling:
%|	A = Gsparse(file)
%| and then you can use it thereafter by typing commands like
%|	y = A * x.
%|
%| in
%|	arg1	char | cell | sparse	sparse matrix (usual case)
%|					or cell array of sparse() arguments
%|					or filename of an aspire .wtf
%| options
%|	mask	[idim]			logical support mask
%|	idim	[1 ndim_in]		input dimensions: nx,ny,nz etc.
%|	odim	[1 ndim_out]		output dimensions: nb,na etc.
%|	chat				verbosity
%|
%| out
%|	ob	[nd np]		nd = prod(odim), np = sum(mask(:))
%|				so it is already "masked"
%|
%| Copyright 2005-6-16, Jeff Fessler, University of Michigan

if nargin == 1 && streq(arg1, 'test'), Gsparse_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

% defaults
arg.mask = [];
arg.idim = [];
arg.odim = [];
%arg.class = 'Fatrix';
arg.class = 'fatrix2';
arg.chat = 0;

arg.blocks = {}; % place to store blocks of G for subset algorithms
arg = vararg_pair(arg, varargin);

if ~isempty(arg.mask) && ~islogical(arg.mask)
	error 'mask must be logical'
end


% cell array of arguments to sparse()
if iscell(arg1)
	if length(arg1) >= 3 % i, j, s ...
		arg1{3} = double(arg1{3}); % trick: double values for s
	end
	arg1 = sparse(arg1{:});
end


% if input is an Aspire .wtf file
if ischar(arg1)
	arg.file = arg1;

	if ~isempty(arg.idim) || ~isempty(arg.odim)
		error 'idim / odim should not be given for .wtf'
	end
	[arg.G arg.idim(1) arg.idim(2) arg.odim(1) arg.odim(2)] = ...
		wtf_read(arg.file);

	% default mask from .wtf
	if isempty(arg.mask)
		tmp = full(sum(arg.G) > 0);
		arg.mask = reshape(tmp, arg.idim);
	end

	arg.G = arg.G(:,arg.mask(:));


% if input is a sparse matrix
elseif issparse(arg1)
	arg.G = arg1;

	if isempty(arg.idim)
		if ~isempty(arg.mask)
			arg.idim = size(arg.mask);
		else
			warn 'idim not given for sparse matrix!'
			arg.idim = [size(arg.G,2) 1];
		end
	end

	if isempty(arg.mask)
		arg.mask = true(arg.idim); % default mask is all
	elseif length(arg.idim) ~= ndims(arg.mask) || ...
		any(arg.idim ~= size(arg.mask))
		disp(arg), error 'bad mask size'
	end

	if isempty(arg.odim)
		warn 'odim not given for sparse matrix!'
		arg.odim = [size(arg.G,1) 1];
	elseif prod(arg.odim) ~= size(arg.G,1)
		error 'bad row dimension'
	end

	if sum(arg.mask(:)) ~= size(arg.G,2)
		if size(arg.G,2) == numel(arg.mask)
			arg.G = arg.G(:, arg.mask(:)); % trick: compact size
		else
			disp(arg), error 'bad G size'
		end
	end

else
	error 'input must be cell or filename or sparse matrix'
end

switch arg.class
case 'Fatrix' % build Fatrix object

	arg.nd = prod(arg.odim);
	arg.np = sum(arg.mask(:));
	ob = Fatrix([arg.nd arg.np], arg, 'caller', mfilename, ...
		'forw', @Gsparse_forw_Fatrix, 'back', @Gsparse_back_Fatrix, ...
		'block_setup', @Gsparse_block_setup, ...
		'mtimes_block', @Gsparse_mtimes_block, ...
		'abs', @Gsparse_abs, 'power', @Gsparse_power);

case 'fatrix2'

	forw = @(arg, x) Gsparse_forw_array(arg, x);
	back = @(arg, y) Gsparse_back_array(arg, y);

	forw_block = @(arg, x, iblock, nblock) ...
		Gsparse_block_forw_array(arg, x, iblock, nblock);
	back_block = @(arg, y, iblock, nblock) ...
		Gsparse_block_back_array(arg, y, iblock, nblock);

	ob = fatrix2('arg', arg, 'mask', arg.mask, 'does_many', 1, ...
		'idim', arg.idim, 'odim', arg.odim, ...
		'block_setup', @Gsparse_block_setup, ...
		'forw', forw, 'back', back, ...
		'forw_block', forw_block, 'back_block', back_block, ...
		'abs', @Gsparse_abs, 'power', @Gsparse_power);

otherwise
	fail 'unknown class'
end


%
% universal versions
%


% Gsparse_abs()
function ob = Gsparse_abs(ob)
ob.arg.G = abs(ob.arg.G);
% fix: inefficient to work with both G and its blocks
if ~isempty(ob.nblock)
	for ii=1:ob.nblock
		ob.arg.blocks{ii} = abs(ob.arg.blocks{ii});
	end
end


% Gsparse_power()
function ob = Gsparse_power(ob, sup)
ob.arg.G = ob.arg.G .^ sup;
% fix: inefficient to work with both G and its blocks
if ~isempty(ob.nblock)
	for ii=1:ob.nblock
		ob.arg.blocks{ii} = ob.arg.blocks{ii} .^ sup;
	end
end


% Gsparse_block_setup()
% Pre-construct blocks of sparse matrix so that it need not be done
% for every block access. Doubles memory.
function ob = Gsparse_block_setup(ob)

nb = prod(ob.arg.odim(1:end-1));
na = ob.arg.odim(end);
ob.arg.blocks = cell(ob.nblock,1);
for iblock=1:ob.nblock
	ia = iblock:ob.nblock:na;
	ii = outer_sum(1:nb, (ia-1)*nb);
	ii = ii(:);
	t = ob.arg.G(ii,:); % fix: sparse, but nzmax is too large!
	if 0 % old matlab version had extra zeros
		[i j s] = find(t);
		if iblock == 1 && length(s) ~= nzmax(t)
%			persistent warned
			warn 'squeezing matlab''s excess zeros in sparse matrix'
		end
		t = sparse(i, j, s, length(ii), size(ob.arg.G, 2));
	end
	ob.arg.blocks{iblock} = t;
end


% Gsparse_forw_mult()
% multiply a sparse matrix by a vector that might not be double
function y = Gsparse_forw_mult(S, x)

if isa(x, 'double')
	y = S * x; % [nd *L]
else
	y = S * double(x); % [nd *L]
	if ~issparse(y)
		y = single(y);
	end
end


% Gsparse_back_mult()
function x = Gsparse_back_mult(S, y)

if isa(y, 'double')
	x = (y' * S)'; % [np *L] trick: runs faster this way!
else
	x = (double(y)' * S)'; % [np *L]
	if ~issparse(x)
		x = single(x);
	end
end
%x = single(full(x)); % todo?


%
% fatrix2 versions
%


% Gsparse_forw_array(): y = A * x
function y = Gsparse_forw_array(arg, x)

x = reshapee(x, numel(arg.mask), []); % [*N *L] array to concise column
x = x(arg.mask(:),:); % [np *L]
y = Gsparse_forw_mult(arg.G, x); % [nd *L]
y = reshaper(y, arg.odim); % [(M) *L]


% Gsparse_back_array(): x = A' * y
function x = Gsparse_back_array(arg, y)

y = reshapee(y, prod(arg.odim), []); % [nd *L]
x = Gsparse_back_mult(arg.G, y); % [np *L]
x = embed(x, arg.mask); % [(N) *L]


% Gsparse_block_forw_array()
function y = Gsparse_block_forw_array(arg, x, istart, nblock);

if istart == 1 && nblock == 1
	y = Gsparse_forw_array(arg, x);
return
end

ia = istart:nblock:arg.odim(end); % subset over last dim

x = reshapee(x, numel(arg.mask), []); % [*N *L] array to concise column
x = x(arg.mask(:),:); % [np *L]
y = Gsparse_forw_mult(arg.blocks{istart}, x); % [nd1 *L]
y = reshaper(y, [arg.odim(1:end-1) length(ia)]);


% Gsparse_block_back_array()
function x = Gsparse_block_back_array(arg, y, istart, nblock);

if istart == 1 && nblock == 1
	x = Gsparse_back_array(arg, y);
return
end

ia = istart:nblock:arg.odim(end); % subset over last dim
nd1 = prod(arg.odim(1:end-1)) * length(ia);

y = reshapee(y, nd1, []); % [nd1 *L]
x = Gsparse_back_mult(arg.blocks{istart}, y); % [np *L]
x = embed(x, arg.mask); % [(N) *L]


%
% Fatrix versions
%


% Gsparse_forw_Fatrix(): y = A * x
function y = Gsparse_forw_Fatrix(arg, x)

% if needed, convert array to concise column
flag_array = 0;
if size(x,1) ~= arg.np
	flag_array = 1;
	x = reshapee(x, numel(arg.mask), []); % [*N *L]
	x = x(arg.mask(:),:); % [np *L]
end

y = Gsparse_forw_mult(arg.G, x); % [nd *L]

if flag_array
	y = reshaper(y, arg.odim); % [(M) *L]
end


% Gsparse_back_Fatrix(): x = A' * y
function x = Gsparse_back_Fatrix(arg, y)

flag_array = 0;
if size(y,1) ~= arg.nd
	flag_array = 1;
	y = reshapee(y, arg.nd, []); % [nd *L]
end

x = Gsparse_back_mult(arg.G, y); % [np *L]

if flag_array
	x = embed(x, arg.mask); % [(N) *L]
end



% Gsparse_mtimes_block()
% note: this is somewhat inefficient, but it is mostly for testing anyway.
function y = Gsparse_mtimes_block(arg, is_transpose, x, istart, nblock)

if is_transpose
	y = Gsparse_mtimes_back(arg, x, istart, nblock);
else
	y = Gsparse_mtimes_forw(arg, x, istart, nblock);
end

% old slow way
%nb = prod(arg.odim(1:end-1));
%ii = outer_sum(1:nb, (ia-1)*nb);
% trick: just reuse almost everything in arg
%arg.G = arg.G(ii(:),:);


% Gsparse_mtimes_forw()
function y = Gsparse_mtimes_forw(arg, x, istart, nblock);

ia = istart:nblock:arg.odim(end); % subset over last dim

% if needed, convert array to concise column
flag_array = 0;
if size(x,1) ~= arg.np
	flag_array = 1;
	x = reshapee(x, numel(arg.mask), []); % [*N *L]
	x = x(arg.mask(:),:); % [np *L]
end

y = Gsparse_forw_mult(arg.blocks{istart}, x); % [nd? *L]

if flag_array
	y = reshaper(y, [arg.odim(1:end-1) length(ia)]);
end


% Gsparse_mtimes_back()
function x = Gsparse_mtimes_back(arg, y, istart, nblock);

ia = istart:nblock:arg.odim(end); % subset over last dim
nd1 = arg.nd * length(ia) / arg.odim(end);

flag_array = 0;
if size(y,1) ~= nd1
	flag_array = 1;
	y = reshapee(y, nd1, []); % [nd1 *L]
end

x = Gsparse_back_mult(arg.blocks{istart}, y); % [np *L]

if flag_array
	x = embed(x, arg.mask); % [(N) *L]
end
