%| argmin(x)
%| argmin(f, x)
%| akin to julia argmin
function out = argmin(f, x)
if isa(f, 'function_handle')
   [~, out] = min(f(x));
else
   [~, out] = min(f);
end
