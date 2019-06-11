 function ob = block_setup(ob, varargin)
%function ob = block_setup(ob, varargin)

if ~isempty(ob.handle_block_setup)
	ob = ob.handle_block_setup(ob, varargin{:});
else
	error(['The Fatrix of type ' ob.caller ' has no handle_block_setup'])
end


% this must be set by the ob.handle_block_setup
%ob.nblock = nblock;

if isempty(ob.nblock) || ob.nblock < 1
	error 'Need nblock >= 1'
end

if isempty(ob.handle_mtimes_block)
	error(['The Fatrix of type ' ob.caller ' has no handle_mtimes_block'])
end

if isempty(ob.handle_blockify_data)
	error(['The Fatrix of type ' ob.caller ' has no handle_blockify_data'])
end
