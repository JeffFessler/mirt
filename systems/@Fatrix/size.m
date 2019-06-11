 function [dim, varargout] = size(ob, varargin)
%function dim = size(ob, varargin)
% or [dim1 dim2] = size(ob)
% "size" method for this class

dim = ob.dim;

if length(varargin) == 1
	dim = dim(varargin{1});
	return
elseif length(varargin) ~= 0
	error 'bug'
end

if nargout > 1
	varargout{1} = dim(2);
	dim = dim(1);
end
