 function [varargout] = dsingle(varargin)
%function y = dsingle(x)
% truncate a (possibly) double precision value to single precision,
% then convert back to "double" using double6:
% y = double6(single(x));

for ii=1:nargin
	varargout{ii} = dsingle1(varargin{ii});
end

function y = dsingle1(x)
y = double6(single(x));
