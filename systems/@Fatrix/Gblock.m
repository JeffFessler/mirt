 function ob = Gblock(ob, nblock, varargin)
%function ob = Gblock(ob, nblock, varargin)

if isempty(ob.handle_mtimes_block) && ob.nblock > 1
	error(['The Fatrix of type ' ob.caller ' has no mtimes_block()'])
end

if ~isfield(ob.arg, 'odim')
	if nblock > 1
		error(['A block object must have arg.odim'])
	else
%		warn 'assuming odim = [nd 1]'
		ob.arg.odim = [ob.dim(1) 1];
	end
end

if nblock > ob.arg.odim(end)
	warn('nblock=%d with dim=[%s]?', nblock, mat2str(ob.arg.odim))
end

ob.nblock = nblock;
ob.caller = sprintf('F:Gblock(%s)', ob.caller);

if ~isempty(ob.handle_block_setup)
	ob = ob.handle_block_setup(ob);
end
