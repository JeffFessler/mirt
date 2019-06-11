 function [str, d1, d2] = poly_string(order, varargin)
%function [str, d1, d2] = poly_string(order, [options])
%|
%| make a string with all the terms in a 2d polynomial
%| up to given order, e.g. [1+0*x x y ...]
%| also return partial derivatives w.r.t. x and y
%|
%| in
%|	order
%|
%| options
%|	maxdegree	maximum of sum of degrees (default: 2*order)
%|	dc		set to 1 to include dc (constant) term (default: 1)
%|
%| out
%|	str		to evaluate polynomial, use:
%|			fun = str2func(str); out = fun(x, y)
%|	d1,d2		likewise, for 1st partial derivatives thereof
%|
%| Copyright 2005-6-18, Jeff Fessler, The University of Michigan

warn 'superceded by ir_poly2_fun'

if nargin < 1, ir_usage, end
if nargin == 1 && streq(order, 'test'), poly_string_test, return, end

arg.maxdegree = [];
arg.dc = true;
arg = vararg_pair(arg, varargin);

if isempty(arg.maxdegree), arg.maxdegree = 2 * order; end

str = '';
d1 = '';
d2 = '';
for i2=0:order
	for i1=0:order
		if i1+i2 == 0 && ~arg.dc
			continue
		end
		if i1 + i2 > arg.maxdegree
			continue
		end

		str = [str sprintf(' (x.^%d).*(y.^%d)', i1, i2)];
		if i1 == 0
			s1 = '0*x';
		else
			s1 = sprintf('%d*x.^%d', i1, i1-1);
		end
		if i2 == 0
			s2 = '0*y';
		else
			s2 = sprintf('%d*y.^%d', i2, i2-1);
		end
		d1 = [d1 sprintf(' (%s).*(y.^%d)', s1, i2)];
		d2 = [d2 sprintf(' (x.^%d).*(%s)', i1, s2)];
	end
end
str = ['[ ' str ' ]'];


% poly_string_test()
function poly_string_test
[str d1 d2] = poly_string(1, 'dc', 0)
[str d1 d2] = poly_string(2, 'maxdegree', 3)
