 function dot = ir_dot_double(a, b)
%function dot = ir_dot_double(a, b)
%| compute dot product between vectors a and b in double precision
%| dot = sum(a(:) .* b(:), 'double'); % double accumulate
%| user must apply conj() to a or b before calling for complex case

if nargin < 2, ir_usage, end

dot = sum(a(:) .* b(:), 'double'); % double accumulate
