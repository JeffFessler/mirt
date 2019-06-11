  function ob = block_diag(varargin)
%|function ob = block_diag(block1, ..., blockM)
%|
%| Construct fatrix2 object from fatrix2 objects
%|	block_diag(A_1, A_2, ..., A_M)
%|
%| in
%|	blocks{:}		fatrix2 blocks
%|
%| options
%|
%| out
%|	ob [nd np]	nd = sum_m nrow(A_M), np = sum_m ncol(A_m)
%|
%| Copyright 2012-09-09, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

% defaults
arg = struct;
% options
% arg = vararg_pair(arg, varargin);

% initialize

ob = fatrix2_block_diag(varargin, arg);


% fatrix2_block_diag()
function ob = fatrix2_block_diag(blocks, arg)

if ~iscell(blocks)
	fail('blocks must be cell array, not "%s"', class(blocks))
end

arg.blocks = blocks;
MM = numel(blocks);
arg.dims = zeros(MM, 2);
for ii=1:MM
	arg.dims(ii,:) = size(blocks{ii});
end
% start/end indices for selecting parts of x and y
arg.istart = cumsum([1; arg.dims(1:end-1,1)]);
arg.iend = arg.istart + arg.dims(:,1) - 1;
arg.jstart = cumsum([1; arg.dims(1:end-1,2)]);
arg.jend = arg.jstart + arg.dims(:,2) - 1;

arg.dim = sum(arg.dims, 1);

% if all the dimensions match
[idims_same idims] = fatrix2_same_idims(blocks);
[odims_same odims] = fatrix2_same_odims(blocks);
if idims_same && odims_same
	arg.idim_cat = 1 + numel(blocks{1}.idim);
	arg.odim_cat = 1 + numel(blocks{1}.odim);
	for ii=1:MM
		imask{ii} = blocks{ii}.imask;
		omask{ii} = blocks{ii}.omask;
	end
	imask = cat(arg.idim_cat, imask{:});
	omask = cat(arg.odim_cat, omask{:});
else
	arg.idim_cat = 0;
	arg.odim_cat = 0;
	fail todo
end

% build object
ob = fatrix2('arg', arg, ...
	'idim', [blocks{1}.idim MM], 'imask', imask, ...
	'odim', [blocks{1}.odim MM], 'omask', omask, ...
	'caller', 'fatrix2_block_diag', ...
	'forw', @fatrix2_block_diag_forw, ...
	'back', @fatrix2_block_diag_back);
%	'free', @fatrix2_block_free, 'gram', @fatrix2_block_diag_gram, ...
%	'mtimes_block', @fatrix2_block_diag_mtimes_block);


% fatrix2_block_diag_forw(): y = A * x
% full forward projection
function y = fatrix2_block_diag_forw(arg, x)

MM = numel(arg.blocks);
if arg.odim_cat
	y = cell(MM, 1);
	for ii=1:MM
		tmp = stackpick(x, ii);
		y{ii} = arg.blocks{ii} * tmp;
	end
	y = cat(arg.odim_cat, y{:});
else
	fail 'todo'
	y = [];
	for ii=1:MM
		t = arg.blocks{ii} * x([arg.jstart(ii):arg.jend(ii)], :);
		y = [y; t];
	end
end


% fatrix2_block_diag_back(): x = A' * y
% full backprojection
function x = fatrix2_block_diag_back(arg, y)

MM = numel(arg.blocks);
if arg.idim_cat
	x = cell(MM, 1);
	for ii=1:MM
		tmp = stackpick(y, ii);
		x{ii} = arg.blocks{ii}' * tmp;
	end
	x = cat(arg.idim_cat, x{:});
else
	fail 'todo'
	x = [];
	for ii=1:MM
		t = arg.blocks{ii}' * y([arg.istart(ii):arg.iend(ii)], :);
		x = [x; t];
	end
end


%{
% fatrix2_block_diag_mtimes_block()
% caution: it is quite unclear whether these are useful!
function y = fatrix2_block_diag_mtimes_block(arg, is_transpose, x, istart, nblock)
if is_transpose
	y = fatrix2_block_diag_block_back(arg, x, istart, nblock);
else
	y = fatrix2_block_diag_block_forw(arg, x, istart, nblock);
end


% fatrix2_block_diag_block_forw()
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
function x = fatrix2_block_diag_block_back(arg, y, istart, nblock)

if nrow(y) ~= arg.dim(1), error 'bad y size', end
x = [];
for ii=istart:nblock:length(arg.blocks)
	t = arg.blocks{ii}' * y([arg.istart(ii):arg.iend(ii)], :);
	x = [x; t];
end
%}
