 function out = equivs(var1, var2, varargin)
%function out = equivs(var1, var2, command)
%|
%| verify that var1 and var2 are equivalent to within single precision accuracy
%| if not, print error message.  an alternative to isequal().
%| See also: jf_equal
%|
%| option
%|	'thresh'	threshold (default: 1e-6)
%|	'format'	format for displaying the two variables (default: '')
%|	'fail'	0|1	if 1 (default) then fail if not equivalent, else warn
%|
%| out
%|	out		0|1
%|
%| Copyright 2007, Jeff Fessler, University of Michigan

if nargin == 1 && streq(var1, 'test'), equivs_test, return, end
if nargin < 2, ir_usage, end

arg.thresh = 1e-6;
arg.format = '';
arg.fail = true;
arg = vararg_pair(arg, varargin);

if isempty(var1) && isempty(var2)
	ok = true;

elseif ~isequal(size(var1), size(var2))
	printm([': size(%s) = %s'], inputname(1), mat2str(size(var1)))
        printm([': size(%s) = %s'], inputname(2), mat2str(size(var2)))
	equivs_show(var1, var2, inputname(1), inputname(2), arg.format)
	if arg.fail
		fail('incompatible dimensions')
	else
		warn('incompatible dimensions')
		return
	end

else
	var1 = var1(:);
	var2 = var2(:);
	if any(isnan(var1)) || any(isnan(var2))
		if arg.fail
			fail('nan')
		else
			warn('nan')
			ok = false;
			return
		end
	end

	norm = (max(abs(var1)) + max(abs(var2))) / 2;
	if ~norm
		ok = true; % both zero!
	else
		err = max(abs(var1-var2)) / norm;
		ok = err < arg.thresh;
	end
end

if nargout
	out = ok;
end

if ok
	return
end

[name line] = caller_name;
if isempty(name)
	str = '';
else
	str = sprintf('%s %d:', name, line);
end

name1 = inputname(1);
name2 = inputname(2);
minmax(var1, ['equivs ' name1 ':'])
minmax(var2, ['equivs ' name2 ':'])
diff = var1 - var2;
if ~isreal(diff)
	diff = abs(diff);
end
minmax(diff)
printm([str ' normalized difference of %g between "%s" "%s"'], ...
	err, name1, name2);
if arg.fail
	fail('not equal to within single precision, thresh=%g', arg.thresh)
else
	warn('not equal to within single precision, thresh=%g', arg.thresh)
end


function equivs_show(var1, var2, name1, name2, format)
if isempty(format), return, end
disp(name1)
disp(num2str(var1, format))
disp(name2)
disp(num2str(var2, format))


function equivs_test
rng(0)
x = randn(1000,200);
y = dsingle(x);
equivs(x,y)

passed = 0;
try
	y = x + 2e-6 * max(x(:));
	equivs(x,y)
	passed = 1;
catch
end
if passed, error 'this should have failed!', end
