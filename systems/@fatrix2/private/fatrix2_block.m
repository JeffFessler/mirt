 function ob = fatrix2_block(blocks, varargin)
%function ob = fatrix2_block(blocks, options)
%|
%| Construct fatrix2_block object, a meta-object composed of fatrix2 blocks,
%| such as block_diag(A_1, A_2, ..., A_M)
%| See fatrix2_block_test.m for example usage.
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

fail 'under construction'

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(blocks, 'test')
	run_mfile_local fatrix2_block_test
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
	ob = fatrix2_block_col(blocks, arg);
	fail 'used fatrix2_vertcat instead'
case 'diag'
	ob = fatrix2_block_diag(blocks, arg);
case 'kron'
	ob = fatrix2_block_kron(blocks, arg);
case 'row'
	ob = fatrix2_block_row(blocks, arg);
case 'sum'
	ob = fatrix2_block_sum(blocks);
otherwise
	fail('unknown block type "%s"', arg.type)
end


% fatrix2_block_row()
% trick: just reuse col via transpose!
%
function ob = fatrix2_block_row(blocks, arg)

for ii=1:length(blocks)
	blocks{ii} = blocks{ii}'; % trick: transpose
end

ob = fatrix2_block(blocks, 'type', 'col')'; % trick: transpose


% fatrix2_block_diag()
%
function ob = fatrix2_block_diag(blocks, arg)

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

% build object
ob = fatrix2(arg.dim, arg, 'caller', 'fatrix2_block(diag)', ...
	'forw', @fatrix2_block_diag_forw, 'back', @fatrix2_block_diag_back, ...
	'free', @fatrix2_block_free, 'gram', @fatrix2_block_diag_gram, ...
	'mtimes_block', @fatrix2_block_diag_mtimes_block);


% fatrix2_block_diag_forw(): y = A * x
%
function y = fatrix2_block_diag_forw(arg, x)

if nrow(x) ~= arg.dim(2)
	error('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
end
y = [];
for ii=1:length(arg.blocks)
	t = arg.blocks{ii} * x([arg.jstart(ii):arg.jend(ii)], :);
	y = [y; t];
end


% fatrix2_block_diag_back(): x = A' * y
% full backprojection
%
function x = fatrix2_block_diag_back(arg, y)

if nrow(y) ~= arg.dim(1), error 'bad y size', end
x = [];
for ii=1:length(arg.blocks)
	t = arg.blocks{ii}' * y([arg.istart(ii):arg.iend(ii)], :);
	x = [x; t];
end


% fatrix2_block_diag_mtimes_block()
% caution: it is quite unclear whether these are useful!
%
function y = fatrix2_block_diag_mtimes_block(arg, is_transpose, x, istart, nblock)
if is_transpose
	y = fatrix2_block_diag_block_back(arg, x, istart, nblock);
else
	y = fatrix2_block_diag_block_forw(arg, x, istart, nblock);
end


% fatrix2_block_diag_block_forw()
%
function y = fatrix2_block_diag_block_forw(arg, x, istart, nblock)

if nrow(x) ~= arg.dim(2)
	error('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
end
y = [];
for ii=istart:nblock:length(arg.blocks)
	t = arg.blocks{ii} * x([arg.jstart(ii):arg.jend(ii)], :);
	y = [y; t];
end


% fatrix2_block_diag_block_back()
%
function x = fatrix2_block_diag_block_back(arg, y, istart, nblock)

if nrow(y) ~= arg.dim(1), error 'bad y size', end
x = [];
for ii=istart:nblock:length(arg.blocks)
	t = arg.blocks{ii}' * y([arg.istart(ii):arg.iend(ii)], :);
	x = [x; t];
end


% fatrix2_block_diag_gram()
%
function [T, reuse] = fatrix2_block_diag_gram(ob, W, reuse, varargin)

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
T = fatrix2_block(T, 'type', 'diag');


% fatrix2_block_kron()
%
function ob = fatrix2_block_kron(blocks, arg)

if isempty(arg.Mkron), error 'Mkron required', end
if length(blocks) ~= 1, error 'kron expects exactly one block', end

arg.dims = repmat(size(blocks{1}), [arg.Mkron 1]);
% start/end indices for selecting parts of x and y
arg.istart = cumsum([1; arg.dims(1:end-1,1)]);
arg.iend = arg.istart + arg.dims(:,1) - 1;
arg.jstart = cumsum([1; arg.dims(1:end-1,2)]);
arg.jend = arg.jstart + arg.dims(:,2) - 1;
arg.dim = sum(arg.dims,1);

% build object
if isa(blocks{1}, 'fatrix2')
	tmp = ['F:' blocks{1}.caller];
else
	tmp = class(blocks{1});
end
tmp = sprintf('fatrix2_block(kron, %s)', tmp);
ob = fatrix2(arg.dim, arg, 'caller', tmp, ...
	'forw', @fatrix2_block_kron_forw, 'back', @fatrix2_block_kron_back, ...
	'fatrix2_block, @fatrix2_block_kron_block_setup, ...
	'mtimes_block', @fatrix2_block_kron_mtimes_block, ...
	'free', @fatrix2_block_free);


% fatrix2_block_kron_forw(): y = A * x
%
function y = fatrix2_block_kron_forw(arg, x)

if nrow(x) ~= arg.dim(2)
	fail('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
end
y = [];
for ii=1:arg.Mkron
	t = arg.blocks{1} * x([arg.jstart(ii):arg.jend(ii)], :);
	y = [y; t];
end


% fatrix2_block_kron_back(): x = A' * y
% full backprojection
%
function x = fatrix2_block_kron_back(arg, y)

if nrow(y) ~= arg.dim(1), error 'bad y size', end
x = [];
for ii=1:arg.Mkron
	t = arg.blocks{1}' * y([arg.istart(ii):arg.iend(ii)], :);
	x = [x; t];
end


% fatrix2_block_kron_block_setup()
% apply Gblock to the base block, to prepare it for later
%
function ob = fatrix2_block_kron_block_setup(ob, varargin)
ob.arg.blocks = {Gblock(ob.arg.blocks{1}, ob.nblock)};
% todo: probably need to set up some other internal variables here
% for use inside fatrix2_block_kron_mtimes_block()


% fatrix2_block_kron_mtimes_block(): y = A{ib} * x
%
function y = fatrix2_block_kron_mtimes_block(arg, is_transpose, x, iblock, nblock)

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
% OLDfatrix2_block_free()
%
function OLDfatrix2_block_free(arg)
if arg.chat
	printm 'freeing fatrix2_block object static memory'
end
for ii=1:length(arg.blocks)
	try
		free(blocks{ii})
	catch
	end
end


%
% fatrix2_block_update()
%
function out = fatrix2_block_update(ob, varargin)
% todo: figure out how to update, e.g., new_zmap(s), ...
