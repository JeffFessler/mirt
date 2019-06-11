 function cells = blockify_data(ob, data, varargin)
%function cells = blockify_data(ob, data, varargin)

if isempty(ob.handle_blockify_data)
	fail('The %s of type %s has no handle_blockify_data', ...
		class(ob), ob.caller)
end
cells = ob.handle_blockify_data(ob, data, varargin{:});
