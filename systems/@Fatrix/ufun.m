 function out = ufun(ob, varargin)
%function out = ufun(ob, varargin)

if ~isempty(ob.handle_ufun)
	out = ob.handle_ufun(ob, varargin{:});
else
	error(['The Fatrix of type ' ob.caller ' has no handle_ufun'])
end
