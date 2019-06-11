 function ob = fatrix2_horzcat(blocks, varargin)
%function ob = fatrix2_horzcat(blocks, [options])
%|
%| horzcat: B = [A1, A2, ...]
%|
%| in
%|	blocks	{cell}	cell array of the (fatrix2) blocks
%|
%| option
%|	'dim_cat'	dimension along which to concatenate idims
%|			all other idims must match
%|			default is ndim+1 if *all* dimensions match
%|			and imask is not full 1D.
%|			otherwise default is 1
%|
%| Copyright 2011-10-02, Jeff Fessler, University of Michigan

arg.dim_cat = []; % default is to infer from idims
arg = vararg_pair(arg, varargin);

MM = numel(blocks);

if 1 % trick: just use transpose of _vertcat
	for ii=1:MM
		blocks{ii} = blocks{ii}'; % trick: transpose
	end

	ob = fatrix2_vertcat(blocks, 'dim_cat', arg.dim_cat)'; % trick: transpose
return
end


% NOT USED BELOW HERE!!!

arg.in_1d = false;

fatrix2_horzcat_check_all_fatrix2(blocks)
fatrix2_horzcat_check_odim(blocks)

[idims_same idims] = fatrix2_same_idims(blocks);
if ~idims_same % columnize
	warn 'mismatched idims so reverting to single column input'
	arg.in_1d = true;
	for mm=1:MM
		idims{mm} = prod(idims{mm});
	end
end

arg.blocks = blocks;

if isempty(arg.dim_cat)
	b1 = blocks{1};
	if numel(b1.odim) == 1 && size(b1,2) == b1.odim % full 1d input
		arg.dim_cat = 1;
	elseif idims_same
		arg.dim_cat = 'next';
	else
		arg.dim_cat = 1;
	end
end


% build object
ob = fatrix2(arg.dim, arg, 'caller', 'fatrix2_horzcat', ...
	'forw', @fatrix2_block_diag_forw, 'back', @fatrix2_block_diag_back, ...
	'free', @fatrix2_block_free, 'gram', @fatrix2_block_diag_gram, ...
	'mtimes_block', @fatrix2_block_diag_mtimes_block);


% fatrix2_horzcat_check_all_fatrix2()
function fatrix2_horzcat_check_all_fatrix2(blocks)
if ~iscell(blocks)
	fail('blocks must be cell array, not "%s"', class(blocks))
end
for mm=1:numel(blocks);
	if ~isa(blocks{mm}, 'fatrix2')
		fail('horzcat requires all fatrix2; use Gmatrix if needed')
	end
end


% fatrix2_horzcat_check_odim()
% check compatibility of output dimensions
function fatrix2_horzcat_check_odim(blocks)
b1 = blocks{1};
for mm=2:numel(blocks);
	bm = blocks{mm};
	if size(bm,1) ~= size(b1,1)
		fail('all blocks must have same #rows')
	end
	if ~isequal(bm.odim, b1.odim)
		fail('all blocks must have same odim')
	end
	if ~isequal(bm.omask, b1.omask)
		fail('all blocks must have same omask')
	end
end


% fatrix2_horzcat_abs(): abs(ob)
function ob = fatrix2_horzcat_abs(ob)

MM = numel(ob.arg.blocks);
for mm=1:MM
	ob.arg.blocks{mm} = abs(ob.arg.blocks{mm});
end


% fatrix2_horzcat_power(): ob.^sup
function ob = fatrix2_horzcat_power(ob, sup)

MM = numel(ob.arg.blocks);
for mm=1:MM
	ob.arg.blocks{mm} = (ob.arg.blocks{mm}) .^ sup;
end
