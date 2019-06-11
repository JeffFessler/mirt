 function h = ir_str2func(string)
%function h = ir_str2func(string)
%|
%| wrapper around str2func() for octave / matlab compatability
%|
%| Background: matlab is phasing out inline
%| http://www.mathworks.com/help/matlab/ref/inline.html
%|
%| So given a string like 'x^2' the way to make it function is
%| str2func('@(x) x^2')
%|
%| Unfortunately, octave 3.8.2 does not support that syntax,
%| so this is the work around.

if nargin < 1, ir_usage, end
if streq(string, 'test') ir_str2func_test, return, end

if ir_is_octave
	h = eval(string);
else
	h = str2func(string);
end

% ir_str2func_test()
% test routine (that is longer than the actual function)
function ir_str2func_test
s = '@(x) x.^2';
f = ir_str2func(s);
g = @(x) x.^2;
jf_equal(f, g)
x = 1:5;
y = f(x);
jf_equal(y, x.^2)
