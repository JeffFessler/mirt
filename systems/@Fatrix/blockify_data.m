 function cells = blockify_data(ob, data, varargin)
%function cells = blockify_data(ob, data, varargin)

if isempty(ob.handle_blockify_data)
	error(['The Fatrix of type ' ob.caller ' has no handle_blockify_data'])
end
cells = ob.handle_blockify_data(ob, data, varargin{:});
