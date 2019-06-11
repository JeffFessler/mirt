function ob = fatrix2_block_sum(blocks)
%| ob = fatrix2_block_sum(blocks)
%| called for ob1 + ob2 + ...

% find at least one fatrix2 in the collection
ib = [];
for ii=1:numel(blocks)
	if isa(blocks{ii}, 'fatrix2')
		ib = ii;
		break
	end
end
if isempty(ib), fail 'no fatrix2? bug', end

bb = blocks{ib}; % reference block

arg.blocks = blocks;

arg.same_idims = true; % retain idim imask if all compatible
arg.same_odims = true; % retain odim omask if all compatible

for ii=1:numel(blocks)
	bi = blocks{ii};
	if ~isequal(size(bi), size(bb))
		error 'blocks must have same size for "sum"'
	end
	if isa(bi, 'fatrix2')
%		try
		bi_imask = fatrix2_imask_array(bi);
		bb_imask = fatrix2_imask_array(bb);
%		catch
%		keyboard
%		end

		if ~isequal(bi_imask, bb_imask)
			if ~isequal(bi_imask(:), bb_imask(:))
				error 'blocks must have same imask(:) for "sum"'
			else
				if arg.same_idims
					warn 'different imask sizes so columnizing'
				end
				arg.same_idims = false;
			end
		end

		if ~isequal(bi.idim, bb.idim)
			if ~isequal(prod(bi.idim), prod(bb.idim))
				error 'blocks must have same prod(idim) for "sum"'
			else
				if arg.same_idims
					warn 'different prod(idim) so columnizing'
				end
				arg.same_idims = false;
			end
		end

		bi_omask = fatrix2_omask_array(bi);
		bb_omask = fatrix2_omask_array(bb);

		if ~isequal(bi_omask, bb_omask)
			if ~isequal(bi_omask(:), bb_omask(:))
				error 'blocks must have same omask(:) for "sum"'
			else
				if arg.same_odims
					warn 'different omask sizes so columnizing'
				end
				arg.same_odims = false;
			end
		end

		if ~isequal(bi.odim, bb.odim)
			if ~isequal(prod(bi.odim), prod(bb.odim))
				error 'blocks must have same prod(odim) for "sum"'
			else
				if arg.same_odims
					warn 'different prod(odim) so columnizing'
				end
				arg.same_odims = false;
			end
		end

	else % force it to be a fatrix!
		arg.blocks{ii} = Gmatrix(bi, ...
			'idim', bb.idim, 'imask', bb.imask, ...
			'odim', bb.odim, 'omask', bb.omask);
	end
end

% build object
if arg.same_idims
	iargs = {'idim', bb.idim, 'imask', bb.imask};
else
	iargs = {'idim', prod(bb.idim), 'imask', bb.imask(:)};
end

if arg.same_odims
	oargs = {'odim', bb.odim, 'omask', bb.omask};
else
	oargs = {'odim', prod(bb.odim), 'omask', bb.omask(:)};
end

ob = fatrix2('caller', 'fatrix2_block(sum)', ...
	'arg', arg, iargs{:}, oargs{:}, ...
	'forw', @fatrix2_block_sum_forw, ...
	'back', @fatrix2_block_sum_back, ...
	'sparse', @fatrix2_block_sum_sparse);
%	'free', @fatrix2_block_free);


% fatrix2_block_sum_forw()
function y = fatrix2_block_sum_forw(arg, x)

for ii=1:numel(arg.blocks)
	bi = arg.blocks{ii};
	tmp = x;
	if ~arg.same_idims
		tmp = reshape(tmp, [bi.idim 1]);
	end
	tmp = fatrix2_do_forw(bi, tmp);
%	tmp = bi.handle_forw(bi.arg, tmp);
%	if bi.scale ~= 1
%		tmp = bi.scale * tmp;
%	end
	if ~arg.same_odims
		tmp = tmp(:);
	end
	if ii == 1
		y = tmp;
	else
		y = y + tmp;
	end
end


% fatrix2_block_sum_back()
function x = fatrix2_block_sum_back(arg, y)

for ii=1:numel(arg.blocks)
	bi = arg.blocks{ii};
	tmp = y;
	if ~arg.same_odims
		tmp = reshape(tmp, [bi.odim 1]);
	end
	tmp = fatrix2_do_back(bi, tmp);
%	tmp = bi.handle_back(bi.arg, tmp);
%	if bi.scale ~= 1
%		tmp = conj(bi.scale) * tmp;
%	end
	if ~arg.same_idims
		tmp = tmp(:);
	end
	if ii == 1
		x = tmp;
	else
		x = x + tmp;
	end
end


% fatrix2_block_sum_sparse()
function sp = fatrix2_block_sum_sparse(arg)

arg = arg.arg; % not sure why

sp = sparse(arg.blocks{1});
for ii=2:numel(arg.blocks)
	sp = sp + sparse(arg.blocks{ii});
end
