  function y = ir_reclass(x, varargin)
%|function y = ir_reclass(x, newclass)
%|
%| Convert input variable x to specific class.
%| Typically: x = ir_reclass(x, class(z)); to match class of z
%|
%| in
%|	x		data
%|
%| option
%|	classname	e.g. 'single' or 'double'
%|
%| out
%|	y		data with converted class
%|
%| 2012-06-04, Jeff Fessler, Univ. of Michigan

if nargin < 1, ir_usage, end
if streq(x, 'test'), ir_reclass_test, return, end

if nargin < 2 % do nothing
	y = x;
return
end

if numel(varargin) ~= 1
	fail 'only one argument allowed'
end

class1 = varargin{1};
switch class1
case 'double'
	y = double(x);
case 'single'
	y = single(x);
otherwise
	fail('class "%s" not done', class1)
end


% ir_reclass_test()
function ir_reclass_test

x = single(7.7);
y = ir_reclass(x, 'double');
equivs(y, x);
y = ir_reclass(x, 'single');
jf_equal(y, x);

x = double(7.7);
y = ir_reclass(x, 'single');
equivs(y, x);
y = ir_reclass(x, 'double');
jf_equal(y, x);
