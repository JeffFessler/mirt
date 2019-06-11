 function y = repeat_slice(x, n)
%function y = repeat_slice(x, n)
%
% kind of a repmat to make 3D from 2D (or 4D from 3D...)

if nargin < 2, ir_usage, end

y = reshape(repmat(x(:), [1 n]), [size(x) n]);
