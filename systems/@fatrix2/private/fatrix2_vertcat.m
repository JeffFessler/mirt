 function ob = fatrix2_vertcat(blocks, varargin)
%function ob = fatrix2_vertcat(blocks, [options])
%| vertcat: B = [A1; A2; ...]
%|
%| in
%|	blocks	{cell}	cell array of the (fatrix2) blocks
%|			all blocks must have same idim and imask
%|
%| option
%|	'dim_cat'	dimension along which to concatenate odims
%|			all other odims must match
%|			default is ndim+1 if *all* dimensions match
%|			and imask is not full 1D.
%|			otherwise default is 1
%| out
%|	ob = [A1; A2; ...]
%|
%| If all input blocks have the same odim and have a handle_forw_block
%| then output could also have a handle_forw_block (todo)

arg.dim_cat = []; % default is to infer from odims
arg = vararg_pair(arg, varargin);
arg.out_1d = false;

MM = numel(blocks);

fatrix2_vertcat_check_all_fatrix2(blocks)
fatrix2_vertcat_check_idim(blocks)

[odims_same, odims] = fatrix2_same_odims(blocks);
if ~odims_same % columnize
	warn 'mismatched odims so reverting to single column output'
	arg.out_1d = true;
	for mm=1:MM
		odims{mm} = prod(odims{mm});
	end
end

arg.blocks = blocks;

if isempty(arg.dim_cat)
	b1 = blocks{1};
	if numel(b1.idim) == 1 && size(b1,2) == b1.idim % full 1d input
		arg.dim_cat = 1;
	elseif odims_same
		arg.dim_cat = 'next';
	else
		arg.dim_cat = 1;
	end
end

if streq(arg.dim_cat, 'next')
	if ~odims_same
		fail 'dim_cat ''next'' works only if all odims same'
	end
	odim = [odims{1} MM];
	arg.dim_cat = numel(odim);
else
	id_check = 1:numel(odims{1});
	id_check(arg.dim_cat) = []; % odim must match except dim_cat
	for mm=1:MM
		if ~isequal(odims{mm}(id_check), odims{1}(id_check))
			fail('all blocks must have same odim (except dim_cat)')
		end
	end

	tmp = cell2mat(odims); % [MM ndim]
	tmp = sum(tmp(:,arg.dim_cat));
	odim = odims{1}; odim(arg.dim_cat) = tmp;
end


if fatrix2_vertcat_omask_all_empty(blocks)
	omask = []; % all empty omasks
else
	tmp = cell(MM,1);
	for mm=1:MM
		tmp{mm} = fatrix2_omask_array(blocks{mm});
		if arg.out_1d
			tmp{mm} = tmp{mm}(:);
		end
	end
	omask = cat(arg.dim_cat, tmp{:});
end

arg.mat2cell_arg = fatrix2_vertcat_mat2cell_arg(odims, arg.dim_cat, odim);

% build object
ob = fatrix2('arg', arg, 'caller', 'fatrix2_vertcat', ...
	'idim', blocks{1}.idim, 'imask', blocks{1}.imask, ...
	'odim', odim, 'omask', omask, ...
	'sparse', @fatrix2_vertcat_sparse, ...
	'abs', @fatrix2_vertcat_abs, 'power', @fatrix2_vertcat_power, ...
	'forw', @fatrix2_vertcat_forw, 'back', @fatrix2_vertcat_back);


% fatrix2_vertcat_check_all_fatrix2()
function fatrix2_vertcat_check_all_fatrix2(blocks)
if ~iscell(blocks)
	fail('blocks must be cell array, not "%s"', class(blocks))
end
for mm=1:numel(blocks);
	if ~isa(blocks{mm}, 'fatrix2')
		fail('vertcat requires all fatrix2; use Gmatrix if needed')
	end
end


% fatrix2_vertcat_check_idim()
% check compatibility of input dimensions
function fatrix2_vertcat_check_idim(blocks)
b1 = blocks{1};
for mm=2:numel(blocks);
	bm = blocks{mm};
	if size(bm,2) ~= size(b1,2)
		fail('all blocks must have same #cols')
	end
	if ~isequal(bm.idim, b1.idim)
		fail('all blocks must have same idim')
	end
	if ~isequal(bm.imask, b1.imask)
		fail('all blocks must have same imask')
	end
end


% fatrix2_vertcat_forw(): y = A * x
function y = fatrix2_vertcat_forw(arg, x)

MM = numel(arg.blocks);
yy = cell(MM,1);
for mm=1:MM
	bi = arg.blocks{mm};
	yy{mm} = fatrix2_do_forw(bi, x);
%	yy{mm} = bi.handle_forw(bi.arg, x);
%	if bi.scale ~= 1
%		yy{mm} = bi.scale * yy{mm};
%	end
	if arg.out_1d
		yy{mm} = yy{mm}(:);
	end
end

y = cat(arg.dim_cat, yy{:});


% fatrix2_vertcat_back(): x = A' * y
function x = fatrix2_vertcat_back(arg, y)

yy = mat2cell(y, arg.mat2cell_arg{:});

MM = length(arg.blocks);
for mm=1:MM
	bi = arg.blocks{mm};
	if arg.out_1d
		yy{mm} = reshape(yy{mm}, [bi.odim 1]);
	end
	tmp = fatrix2_do_back(bi, yy{mm});
%	tmp = bi.handle_back(bi.arg, yy{mm});
%	if bi.scale ~= 1
%		tmp = conj(bi.scale) * tmp;
%	end
	if mm == 1
		x = tmp;
	else
		x = x + tmp;
	end
end


% fatrix2_vertcat_abs(): abs(ob)
function ob = fatrix2_vertcat_abs(ob)

MM = numel(ob.arg.blocks);
for mm=1:MM
	ob.arg.blocks{mm} = abs(ob.arg.blocks{mm});
end


% fatrix2_vertcat_power(): ob.^sup
function ob = fatrix2_vertcat_power(ob, sup)

MM = numel(ob.arg.blocks);
for mm=1:MM
	ob.arg.blocks{mm} = (ob.arg.blocks{mm}) .^ sup;
end


% fatrix2_vertcat_mat2cell_arg()
function mat2cell_arg = fatrix2_vertcat_mat2cell_arg(odims, dim_cat, odim)

MM = numel(odims); % odims is cell(MM,1)

if dim_cat == numel(odims{1}) + 1
	for id = 1:numel(odims{1})
		mat2cell_arg{id} = odim(id);
	end
	mat2cell_arg{dim_cat} = ones(MM,1);
	tmp1 = 1:prod(odim);
	tmp1 = reshape(tmp1, odim);
	tmp2 = mat2cell(tmp1, mat2cell_arg{:});
else
	for mm=1:MM
		odim_m = odims{mm};
		if mm==1
			for id = 1:numel(odim_m)
				if id ~= dim_cat
					mat2cell_arg{id} = odim_m(id);
				end
			end
		end
		mat2cell_arg{dim_cat}(mm) = odim_m(dim_cat);
	end
end


% fatrix2_vertcat_omask_all_empty()
% see if all blocks have empty omask (often this will be the case)
function allempty = fatrix2_vertcat_omask_all_empty(blocks)
allempty = true;
for mm=1:numel(blocks)
	allempty = allempty & isempty(blocks{mm}.omask);
end


% fatrix2_vertcat_sparse(): S = sparse(A)
function S = fatrix2_vertcat_sparse(ob)
arg = ob.arg;
S = sparse([]);
MM = numel(arg.blocks);
S = cell(MM,1);
for mm=1:MM
	S{mm} = sparse(arg.blocks{mm});
end
S = cat(1, S{:});


%{
% fatrix2_vertcat_gram()
% todo: this would be useful if W is empty or (block) diagonal
% because it would save memory, but it is extra complexity
function [T, reuse] = fatrix2_vertcat_gram(ob, W, reuse, varargin)

blocks = ob.arg.blocks;
T = cell(size(blocks));
for mm=1:length(blocks)
	A = blocks{mm};
	if isnumeric(A)
		if isempty(W)
			T{mm} = A' * A;
		else
			warn 'todo: this may not work, need piece of W!'
			T{mm} = A' * W * A;
		end
	else
		if isempty(W)
			T{mm} = build_gram(A, [], reuse, varargin{:});
		else
			if isvar('W.arg.blocks{mm}')
				T{mm} = build_gram(A, W.arg.blocks{mm}, ...
					reuse, varargin{:});
			else
				fail('fatrix2_vertcat_gram needs block diag W')
			end
		end
	end
end
T = fatrix_block_sum(T{:});
%}
