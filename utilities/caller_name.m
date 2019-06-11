  function [name, line] = caller_name(level)
%|function [name, line] = caller_name(level)
%|
%| return name (and line) of calling routine or file (if level=1, the default),
%| or name further up or down the stack by changing level.
%|
%| freemat does not have 'dbstack' so return placeholders: name='?' and line=-1

if nargin < 1
	level = 1;
end

if ~exist('dbstack', 'builtin') % freemat
	name = '?';
	line = -1;
else
	st = dbstack;
	if length(st) < 2+level
		name = '';
		line = 0;
	else
		st = st(2+level);
		name = st.name;
		line = st.line;
	end
end
