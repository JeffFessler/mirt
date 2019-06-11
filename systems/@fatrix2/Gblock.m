 function ob = Gblock(ob, nblock, varargin)
%function ob = Gblock(ob, nblock, varargin)
%|
%| Create "block" fatrix (blocks of rows) for block-based methods like OSEM
%| option:
%|	'idiag'	[idim]	optional input diagonal scaling matrix
%|	'odiag'	[odim]	optional output diagonal scaling matrix

arg.idiag = [];
arg.odiag = [];
arg = vararg_pair(arg, varargin);

if ~isempty(arg.idiag)
	if ~isequal(size(arg.idiag), ob.idim)
		pr size(arg.idiag)
		pr ob.idim
		fail('bad idiag size')
	end
	ob.idiag = arg.idiag; % save optional "input" diagonal scaling matrix
end

if ~isempty(arg.odiag)
	if ~isequal(size(arg.odiag), ob.odim)
		pr size(arg.odiag)
		pr ob.odim
		fail('bad odiag size')
	end
	ob.odiag = arg.odiag; % save optional "output" diagonal scaling matrix
end

if (nblock > 1) && ...
	(isempty(ob.handle_forw_block) || isempty(ob.handle_back_block))
%	isempty(ob.handle_mtimes_block) &
	fail('The %s of type %s is missing forw/back_block()', ...
		class(ob), ob.caller)
end

if size(ob,1) ~= prod(ob.odim)
	fail 'block fatrix2 needs size(ob,1) == prod(odim)'
end

if nblock == 1
	ob.odim = [ob.odim 1];
end

if nblock > ob.odim(end)
	warn('nblock=%d with dim=[%s]?', nblock, mat2str(ob.odim))
end

ob.nblock = nblock;
ob.caller = sprintf('f:Gblock(%s)', ob.caller);

if ~isempty(ob.handle_block_setup)
	ob = ob.handle_block_setup(ob);
end
