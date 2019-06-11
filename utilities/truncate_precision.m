 function y = truncate_precision(x, digits, varargin)
%function y = truncate_precision(x, digits, [option])
%|
%| truncate x to "digits" significant digits of precision
%| option
%|	'type'		floor (default) | round | ceil
%|			floor means "towards 0" and ceil means "away from 0"
%|
%| Copyright 2008-11-05, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), truncate_precision_test, return, end
if nargin < 2, ir_usage, end

arg.type = 'floor';
arg = vararg_pair(arg, varargin);

sgn = sign(x);
x = abs(x);
x(sgn == 0) = 1; % avoid warnings
pow = log10(x);
pow = floor(pow);
x = x ./ 10.^pow;
x = x .* 10.^(digits-1);
switch arg.type
case 'floor'
	x = floor(x);
case 'ceil'
	x = ceil(x);
case 'round'
	x = round(x);
otherwise
	fail('bad type %s', arg.type)
end
x = x ./ 10.^(digits-1);
y = x .* 10.^pow;

y(sgn == 0) = 0;
y = y .* sgn;


% truncate_precision_test()
function truncate_precision_test
y = truncate_precision([1234 0.1234 5 0 -1299], [3 2 2 2 2]);
jf_equal(y, [1230 0.12 5 0 -1200])
y = truncate_precision([-pi pi], 3, 'type', 'ceil');
jf_equal(y, [-3.15 3.15]);
y = truncate_precision([-pi pi], 3, 'type', 'floor');
jf_equal(y, [-3.14 3.14]);
y = truncate_precision(255+1e-6, 3); % todo: not happy with this case
equivs(y, 255)
