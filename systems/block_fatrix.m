  function ob = block_fatrix(blocks, varargin)
%|function ob = block_fatrix(blocks, options)
%|
%| Construct block_fatrix object, a meta-Fatrix composed of Fatrix blocks,
%| such as block_diag(A_1, A_2, ..., A_M)
%| See block_fatrix_test.m for example usage.
%|
%| in
%|	blocks	{cell}	cell array of the blocks
%|
%| options
%|	'type'	char	options:
%|				'diag' for block diagonal (default)
%|				'col' for [A_1; A_2; ...; A_M]
%|				'kron' for kron(eye(Mkron), blocks{1})
%|				'row' for [A1, A2, ..., A_M]
%|				'sum' for A1 + A2 + ... + A_M
%|	'Mkron'	int	required for 'kron' type
%|	'tomo'	0|1	special support for tomography-type objects (default: 0)
%|
%| out
%|	ob [nd np]	diag,row: np = sum_m ncol(A_m)
%|
%| Copyright 05-5-12, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(blocks, 'test')
	run_mfile_local block_fatrix_test
return
end

arg.blocks = blocks;

% defaults
arg.type = 'diag';
arg.chat = 0;
arg.tomo = 0;
arg.Mkron = [];

% options
arg = vararg_pair(arg, varargin);

% initialize
if ~iscell(blocks)
	fail('blocks must be cell array, not "%s"', class(blocks))
end

switch arg.type
case 'col'
	ob = block_fatrix_col(blocks, arg);
case 'diag'
	ob = block_fatrix_diag(blocks, arg);
case 'kron'
	ob = block_fatrix_kron(blocks, arg);
case 'row'
	ob = block_fatrix_row(blocks, arg);
case 'sum'
	ob = block_fatrix_sum(blocks, arg);
otherwise
	fail('unknown block type "%s"', arg.type)
end


%
% block_fatrix_col()
%
function ob = block_fatrix_col(blocks, arg)

MM = length(blocks);
arg.dims = zeros(MM, 2);
for ii=1:MM
	arg.dims(ii,:) = size(blocks{ii});
	if arg.dims(ii,2) ~= arg.dims(1,2)
		error 'all blocks must have same #cols for "col"'
	end
end
% start/end indices for selecting parts of x and y
arg.istart = cumsum([1; arg.dims(1:end-1,1)]);
arg.iend = arg.istart + arg.dims(:,1) - 1;

arg.dim = [sum(arg.dims(:,1)) arg.dims(1,2)];

if arg.tomo % prep for mat2cell later
	for ii=1:MM
		odim = arg.blocks{ii}.arg.odim;
		arg.tomo_ndim = length(odim);
		arg.tomo_dim = arg.tomo_ndim; % insist on last dimension
		if ii==1
			for id = 1:arg.tomo_ndim
				if id ~= arg.tomo_dim
					arg.mat2cell_arg{id} = odim(id);
				end
			end
		end
		arg.mat2cell_arg{arg.tomo_dim}(ii) = odim(arg.tomo_dim);
	end
end

% build Fatrix object
ob = Fatrix(arg.dim, arg, 'caller', 'block_fatrix(col)', ...
	'forw', @block_fatrix_col_forw, 'back', @block_fatrix_col_back, ...
	'free', @block_fatrix_free, 'gram', @block_fatrix_col_gram);


%
% block_fatrix_col_forw(): y = A * x
%
function y = block_fatrix_col_forw(arg, x)

MM = length(arg.blocks);
y = cell(MM,1);
for ii=1:MM
	y{ii} = arg.blocks{ii} * x;
end

if ~arg.tomo || ncol(y{1}) == 1 ...
	|| (arg.tomo && arg.dim(2) == size(x,1)) % [np (L)]
	y = cat(1, y{:});
else % handle 'tomo' case where x is not a single column
	y = cat(arg.tomo_dim, y{:});
end


%
% block_fatrix_col_back(): x = A' * y
%
function x = block_fatrix_col_back(arg, y)

MM = length(arg.blocks);
if arg.tomo
	if arg.dim(1) == size(y,1) % [nd (L)]
		arg.tomo = 0;
	else % [(odim) (L)] split into separate parts, undoing cat(?, y{:})
		tmp = size(y);
		if length(tmp) == arg.tomo_ndim
			yc = mat2cell(y, arg.mat2cell_arg{:});
		elseif length(tmp) > arg.tomo_ndim
			yc = mat2cell(y, arg.mat2cell_arg{:}, ...
				tmp((arg.tomo_ndim+1):end)); % handle (L)
		else
			error 'bug'
		end
	end
end

x = 0;
for ii=1:MM
	if ~arg.tomo % [nd (L)]
		jj = arg.istart(ii):arg.iend(ii);
		yi = y(jj,:);
	else
		yi = yc{ii};
	end
	t = arg.blocks{ii}' * yi;
	x = x + t;
end


%
% block_fatrix_col_gram()
%
function [T, reuse] = block_fatrix_col_gram(ob, W, reuse, varargin)

blocks = ob.arg.blocks;
T = cell(size(blocks));
for ii=1:length(blocks)
	A = blocks{ii};
	if isnumeric(A)
		if isempty(W)
			T{ii} = A' * A;
		else
			warn 'todo: this may not work, need piece of W!'
			T{ii} = A' * W * A;
		end
	else
		if isempty(W)
			T{ii} = build_gram(A, [], reuse, varargin{:});
		else
			if isvar('W.arg.blocks{ii}')
				T{ii} = build_gram(A, W.arg.blocks{ii}, ...
					reuse, varargin{:});
			else
				fail('block_fatrix_col_gram needs block diag W')
			end
		end
	end
end
T = fatrix_plus(T{:});


%
% block_fatrix_row()
% trick: just reuse col via transpose!
%
function ob = block_fatrix_row(blocks, arg)

for ii=1:length(blocks)
	blocks{ii} = blocks{ii}'; % trick: transpose
end

ob = block_fatrix(blocks, 'type', 'col')'; % trick: transpose


%
% block_fatrix_diag()
%
function ob = block_fatrix_diag(blocks, arg)

arg.dims = zeros(length(blocks), 2);
for ii=1:length(blocks)
	arg.dims(ii,:) = size(blocks{ii});
end
% start/end indices for selecting parts of x and y
arg.istart = cumsum([1; arg.dims(1:end-1,1)]);
arg.iend = arg.istart + arg.dims(:,1) - 1;
arg.jstart = cumsum([1; arg.dims(1:end-1,2)]);
arg.jend = arg.jstart + arg.dims(:,2) - 1;

arg.dim = sum(arg.dims, 1);

% build Fatrix object
ob = Fatrix(arg.dim, arg, 'caller', 'block_fatrix(diag)', ...
	'forw', @block_fatrix_diag_forw, 'back', @block_fatrix_diag_back, ...
	'free', @block_fatrix_free, 'gram', @block_fatrix_diag_gram, ...
	'mtimes_block', @block_fatrix_diag_mtimes_block);


%
% block_fatrix_diag_forw(): y = A * x
%
function y = block_fatrix_diag_forw(arg, x)

if nrow(x) ~= arg.dim(2)
	error('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
end
y = [];
for ii=1:length(arg.blocks)
	t = arg.blocks{ii} * x([arg.jstart(ii):arg.jend(ii)], :);
	y = [y; t];
end


%
% block_fatrix_diag_back(): x = A' * y
% full backprojection
%
function x = block_fatrix_diag_back(arg, y)

if nrow(y) ~= arg.dim(1), error 'bad y size', end
x = [];
for ii=1:length(arg.blocks)
	t = arg.blocks{ii}' * y([arg.istart(ii):arg.iend(ii)], :);
	x = [x; t];
end


%
% block_fatrix_diag_mtimes_block()
% caution: it is quite unclear whether these are useful!
%
function y = block_fatrix_diag_mtimes_block(arg, is_transpose, x, istart, nblock)
if is_transpose
	y = block_fatrix_diag_block_back(arg, x, istart, nblock);
else
	y = block_fatrix_diag_block_forw(arg, x, istart, nblock);
end


%
% block_fatrix_diag_block_forw()
%
function y = block_fatrix_diag_block_forw(arg, x, istart, nblock)

if nrow(x) ~= arg.dim(2)
	error('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
end
y = [];
for ii=istart:nblock:length(arg.blocks)
	t = arg.blocks{ii} * x([arg.jstart(ii):arg.jend(ii)], :);
	y = [y; t];
end


%
% block_fatrix_diag_block_back()
%
function x = block_fatrix_diag_block_back(arg, y, istart, nblock)

if nrow(y) ~= arg.dim(1), error 'bad y size', end
x = [];
for ii=istart:nblock:length(arg.blocks)
	t = arg.blocks{ii}' * y([arg.istart(ii):arg.iend(ii)], :);
	x = [x; t];
end


%
% block_fatrix_diag_gram()
%
function [T, reuse] = block_fatrix_diag_gram(ob, W, reuse, varargin)

blocks = ob.arg.blocks;
T = cell(size(blocks));
for ii=1:length(blocks)
	A = blocks{ii};
	if isnumeric(A)
		if isempty(W)
			T{ii} = A' * A;
		else
			T{ii} = A' * W * A;
		end
	else
		T{ii} = build_gram(A, W, reuse, varargin{:});
	end
end
T = block_fatrix(T, 'type', 'diag');



%
% block_fatrix_kron()
%
function ob = block_fatrix_kron(blocks, arg)

if isempty(arg.Mkron), error 'Mkron required', end
if length(blocks) ~= 1, error 'kron expects exactly one block', end

arg.dims = repmat(size(blocks{1}), [arg.Mkron 1]);
% start/end indices for selecting parts of x and y
arg.istart = cumsum([1; arg.dims(1:end-1,1)]);
arg.iend = arg.istart + arg.dims(:,1) - 1;
arg.jstart = cumsum([1; arg.dims(1:end-1,2)]);
arg.jend = arg.jstart + arg.dims(:,2) - 1;
arg.dim = sum(arg.dims,1);

% build Fatrix object
if isa(blocks{1}, 'Fatrix')
	tmp = ['F:' blocks{1}.caller];
else
	tmp = class(blocks{1});
end
tmp = sprintf('block_fatrix(kron, %s)', tmp);
ob = Fatrix(arg.dim, arg, 'caller', tmp, ...
	'forw', @block_fatrix_kron_forw, 'back', @block_fatrix_kron_back, ...
	'block_setup', @block_fatrix_kron_block_setup, ...
	'mtimes_block', @block_fatrix_kron_mtimes_block, ...
	'free', @block_fatrix_free);


%
% block_fatrix_kron_forw(): y = A * x
%
function y = block_fatrix_kron_forw(arg, x)

if nrow(x) ~= arg.dim(2)
	fail('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
end
y = [];
for ii=1:arg.Mkron
	t = arg.blocks{1} * x([arg.jstart(ii):arg.jend(ii)], :);
	y = [y; t];
end


%
% block_fatrix_kron_back(): x = A' * y
% full backprojection
%
function x = block_fatrix_kron_back(arg, y)

if nrow(y) ~= arg.dim(1), error 'bad y size', end
x = [];
for ii=1:arg.Mkron
	t = arg.blocks{1}' * y([arg.istart(ii):arg.iend(ii)], :);
	x = [x; t];
end


%
% block_fatrix_kron_block_setup()
% apply Gblock to the base block, to prepare it for later
%
function ob = block_fatrix_kron_block_setup(ob, varargin)
ob.arg.blocks = {Gblock(ob.arg.blocks{1}, ob.nblock)};
% todo: probably need to set up some other internal variables here
% for use inside block_fatrix_kron_mtimes_block()


%
% block_fatrix_kron_mtimes_block(): y = A{ib} * x
%
function y = block_fatrix_kron_mtimes_block(arg, is_transpose, x, iblock, nblock)

bl1 = arg.blocks{1}; % base block, already put through Gblock
warn 'todo: size check not done'

if ~is_transpose % forw
%	if nrow(x) ~= arg.dim(2)
%		fail('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
%	end
	y = [];
	for ii=1:arg.Mkron
		t = bl1{iblock} * x([arg.jstart(ii):arg.jend(ii)], :);
		y = [y; t];
	end

else % back
	y = x;
%	if nrow(y) ~= arg.dim(1), error 'bad y size', end
	x = [];
	for ii=1:arg.Mkron
		tmp = [arg.istart(ii):arg.iend(ii)]; % i list (if all data)
		% todo: we need a certain subset of that list
		fail('todo: kron subset backprojector not done')
		t = bl1{iblock}' * y(tmp, :);
		x = [x; t];
	end
	y = x;
end


%
% block_fatrix_sum()
%
function ob = block_fatrix_sum(blocks, arg)

dim = size(blocks{1});
for ii=1:length(blocks)
	if ~isequal(size(blocks{ii}), dim)
		error 'blocks must have same size for "sum"'
	end
end

% build Fatrix object
ob = Fatrix(dim, arg, 'caller', 'block_fatrix(sum)', ...
	'forw', @block_fatrix_sum_forw, 'back', @block_fatrix_sum_back, ...
	'free', @block_fatrix_free);


%
% block_fatrix_sum_forw(): y = A * x
%
function y = block_fatrix_sum_forw(arg, x)

y = 0;
for ii=1:length(arg.blocks)
	y = y + arg.blocks{ii} * x;
end


%
% block_fatrix_sum_back(): x = A' * y
% full backprojection
%
function x = block_fatrix_sum_back(arg, y)

x = 0;
for ii=1:length(arg.blocks)
	x = x + arg.blocks{ii}' * y;
end


%
% block_fatrix_free()
%
function block_fatrix_free(arg)
if arg.chat
	printm 'freeing block_fatrix object static memory'
end
for ii=1:length(arg.blocks)
	try
		free(blocks{ii})
	catch
	end
end


%
% block_fatrix_update()
%
function out = block_fatrix_update(ob, varargin)
% todo: figure out how to update, e.g., new_zmap(s), ...
