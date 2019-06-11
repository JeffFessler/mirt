  function [out ok] = ir_struct_find_field(ss, field)
%|function [out ok] = ir_struct_find_field(ss, field)
%|
%| Look recursively through the fields of a structure for a given field.
%| Objects like Fatrix and fatrix2 are converted to structs too for this.
%|
%| Returns ss.(field) or ss.?.(field) or ss.?.?.(field) etc.

if nargin == 1 && streq(ss, 'test'), ir_struct_find_field_test, return, end
if nargin < 2, ir_usage, end

[out ok] = ir_struct_find_field_do(ss, field);
if ~ok
	fail('could not find field "%s" in "%s"', field, inputname(1))
end


function [out ok] = ir_struct_find_field_do(ss, field)

out = [];
ok = false;

if isfield(ss, field)
	out = ss.(field);
	ok = true;
return
end

fnames = fieldnames(ss);
for ii=1:length(fnames)
	fname = fnames{ii};

	s1 = ss.(fname);
	if isa(s1, 'Fatrix') || isa(s1, 'fatrix2')
		s1 = struct(s1); % trick
	end

	if isstruct(s1)
		[out ok] = ir_struct_find_field_do(s1, field);
		if ok
			return
		end
	end
end


% ir_struct_find_field_test()
function ir_struct_find_field_test
ss.a = 1;
ss.b.c = 2;
ss.b.d = 3;
jf_equal(1, ir_struct_find_field(ss, 'a'))
jf_equal(3, ir_struct_find_field(ss, 'd'))
try
	ir_struct_find_field(ss, 'bad')
	fail 'should not get here'
catch
end
