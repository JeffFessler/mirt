  function yn = jf_equal(a, b, varargin)
%|function yn = jf_equal(a, b, varargin)
%|
%| verify that the two arguments are equal.
%| if not, print informative error message.
%| See also: equivs
%| note: to compare objects, use jf_equal(struct(a), struct(b))
%|
%| option
%|	'what'			'warn' 'fail' 'silent'
%|				default: 'silent' if nargout, 'fail' otherwise
%|	'name1'		''	name for input argment a (default: inputname(1))
%|	'name2'		''	name for input argment a (default: inputname(2))
%|	'func2str'	0|1	convert function_handle to char? (default: 1)
%|	'accept_logical_eq_uint8'	(default: false) equate logical with uint8?
%|
%| out
%|	yn	1|0	1 if equal
%|
%| Copyright 2007, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(a, 'test'), jf_equal_test, return, end

arg.what = '';
arg.func2str = true;
arg.name1 = inputname(1);
arg.name2 = inputname(2);
arg.accept_logical_eq_uint8 = false; % equate logical with uint8?
arg = vararg_pair(arg, varargin);

if isempty(arg.what)
	if nargout
		arg.what = 'silent';
	else
		arg.what = 'fail';
	end
end

switch arg.what
case 'silent'
	fun = @do_nothing;
case 'warn'
	fun = @warn;
case 'fail'
	fun = @fail;
otherwise
	fail('bad what = "%s"', arg.what)
end

[name line] = caller_name;
if isempty(name)
	str = '';
else
	str = sprintf('%s %d', name, line);
end

if streq(class(a), 'strum') || streq(class(a), 'fatrix2')
	a = struct(a);
end

if streq(class(b), 'strum') || streq(class(b), 'fatrix2')
	b = struct(b);
end

if islogical(a) && isa(b, 'uint8') && arg.accept_logical_eq_uint8
	a = uint8(a);
end

if islogical(b) && isa(a, 'uint8') && arg.accept_logical_eq_uint8
	b = uint8(b);
end

if (~isnumeric(a) || ~isnumeric(b)) && ~isequal(class(a), class(b))
	if nargout
		yn = false;
	end
	fun([str ': "%s" (%s) and "%s" (%s) unequal class'], ...
		arg.name1, class(a), arg.name2, class(b))
return
end

% trick: matlab calls identical function handles unequal, so convert to char
if arg.func2str && isa(a, 'function_handle')
	a = func2str(a);
	b = func2str(b);
end


if isstruct(a)
	arg_what = arg.what;
	if streq(arg.what, 'fail')
		arg_what = 'warn'; % trick to descend all 
	end
	yn = jf_equal_struct(a, b, arg.name1, arg.name2, arg_what, arg.func2str);
else
%	isequalwithequalnans % could use this instead?
	yn = isequal(a, b);
end

if ~yn && ~streq(arg.what, 'silent')
	jf_not_equal_details(a, b, arg.name1, arg.name2, str, fun)
end
if ~nargout
	clear yn
end


% jf_not_equal_details()
function jf_not_equal_details(a, b, name1, name2, str, fun)

switch class(a)
case 'struct'
	fun([str '"%s" and "%s" unequal struct'], name1, name2)

case 'char'
	fun([str ': "%s" (%s) and "%s" (%s) unequal char'], ...
		name1, a, name2, b)

case 'function_handle'
	a = func2str(a);
	b = func2str(b);
	if ~isequal(a,b)
		fun([str ': "%s" (%s) and "%s" (%s) unequal function_handle'], ...
			name1, a, name2, b)
	end

otherwise

	if islogical(a) && islogical(b)
		a = uint8(a); b = uint8(b); % make "numeric" to compare
	end

	if isnumeric(a) && isnumeric(b) % case {'double', 'single'}
		if max(size(a)) == 1 && max(size(b)) == 1
			fun([str ': "%s" (%g) and "%s" (%g) unequal'], ...
				name1, a, name2, b)
		else
			minmax(a)
			minmax(b)
			if ~isequal(size(a), size(b))
				printm(['size(%s) = %s'], name1, mat2str(size(a)))
				printm(['size(%s) = %s'], name2, mat2str(size(b)))
				error 'dimension mismatch'
			end
			max_percent_diff(a, b)
			fun([str ': "%s" and "%s" unequal'], name1, name2)
		end

	else
		fun([str ': "%s" (%s) and "%s" (%s) unknown class'], ...
			name1, class(a), name2, class(b))
	end
end


% jf_equal_struct()
% see if two structs are equal, descending into fields recursively
function yn = jf_equal_struct(s1, s2, name1, name2, arg_what, arg_func2str)
f1 = fieldnames(s1);
f2 = fieldnames(s2);
yn = ~xor(isempty(f1), isempty(f2));
if isempty(f1) || isempty(f2)
	return
end

try
	s2 = orderfields(s2, s1);
catch
	if streq(arg_what, 'warn')
		warn('%s %s have different fields', name1, name2)
	end
	yn = false;
return
end

names = fieldnames(s1);
for ii=1:length(names)
	name = names{ii};
	f1 = getfield(s1, name);
	f2 = getfield(s2, name);
	yn = yn && jf_equal(f1, f2, 'what', arg_what, ...
			'func2str', arg_func2str, ...
			'name1', [name1 '.' name], ...
			'name2', [name2 '.' name]);
end


% do_nothing()
function do_nothing(varargin)


function jf_equal_test
a = 7;
b = 7;
c = 8;
jf_equal(a,b)
jf_equal(a,7)
%jf_equal(a,c)

if 1 % test struct case
	x.a = 7;
	x.b = '8';
	x.c.d = 3;
	if 1 % equal struct
		y = x;
	else
		y.a = 8;
		y.b = '9';
		y.c.d = 4;
	end
	jf_equal(x, y)
end
