 function y = reale(x, arg2, arg3)
%| Return real part of complex data (with error checking).
%function y = reale(x, arg2, arg3)
%|
%| y = reale(x)
%| y = reale(x, tol)
%| y = reale(x, 'warn', 'message')
%| y = reale(x, 'error')
%| y = reale(x, 'report')
%| y = reale(x, 'prompt')
%| y = reale(x, 'disp')
%|
%| Checks that imaginary part is negligible (warns etc. if not).
%|
%| Copyright Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(x, 'test'), reale_test, return, end

com = 'error';
if isa(x, 'double')
	tol = 1e-13;
else
	tol = 1e-6;
end

if nargin > 1
	if ischar(arg2)
		com = arg2;
	elseif isnumeric(arg2)
		tol = arg2;
	end
end

switch com
case {'disp','prompt', 'report'}
	;
case 'warn'
	onlywarn = 1;
case 'error'
	onlywarn = 0;
otherwise
	fail('bad argument "%s"', com)
end

max_abs_x = max(abs(x(:)));
if max_abs_x == 0
	if any(imag(x(:)) ~= 0)
		fail 'max real 0, but imaginary!'
	else
		y = real(x);
		return
	end
end

frac = max(abs(imag(x(:)))) / max_abs_x;
if streq(com, 'report')
	printm('imaginary part %g%%', frac * 100)
	return
end

if frac > tol
	[cname line] = caller_name;
	t = sprintf('%s(%d): %s: imaginary fraction of %s [class %s] is %g', ...
		cname, line, mfilename, inputname(1), class(x), frac);
	if isvar('arg3')
		t = [t ', ' arg3];
	end
	if streq(com, 'disp')
		disp(t)

	elseif streq(com, 'prompt')
		printm('reale() called for input with imaginary part %g%%', frac * 100)
		printm('reale() called in context where a large imaginary part')
		printm('is likely an *error*.  proceed with caution!')
		t = input('proceed? [y|n]: ', 's');
		if isempty(t) || t(1) ~= 'y'
			printm('ok, aborting is probably wise!')
			error ' '
		end

	elseif onlywarn
		disp(t)
	else
		fail(t)
	end
end
y = real(x);


function reale_test
x = 7 + 1i*eps;
reale(x, 'warn');
reale(x, 'prompt');
%reale(x, 'report'); % check error reporting
%reale(x, eps/100) % check error faulting
