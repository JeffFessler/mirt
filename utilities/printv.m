 function printv(arg)
%function printv(arg)
% print a scalar variable
% todo: add 'off' and 'on'
if nargin < 1, ir_usage, end

printf('"%s" = %g', inputname(1), arg)
