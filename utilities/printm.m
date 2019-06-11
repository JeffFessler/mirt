 function printm(varargin)
%function printm(varargin)
%|
%| like printf except that it puts the mfilename in front of it
%| so that you know where the message originated.

[caller, line] = caller_name;
if ~isempty(line) && line ~= 0
	caller = [caller sprintf(' %d', line)];
	caller = strrep(caller, 'LiveEditorEvaluationHelper', 'LiveE*');
end
if length(varargin)
	tmp = varargin{1};
	if streq(tmp, '\n', 2) % trick for newline to precede caller name
		disp(' ') % blank line
		varargin{1} = tmp(3:end); % strip '\n' at start
	end
	disp([caller ': ' sprintf(varargin{:})])
else
	disp([caller ': '])
end
if ir_is_octave
	fflush(stdout);
end
