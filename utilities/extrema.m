%| extrema(x)
%| akin to julia extrema
function out = extrema(f, x)
if isa(f, 'function_handle')
   out = minmax(f(x));
else
   out = minmax(f);
end
