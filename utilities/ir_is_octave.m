  function out = ir_is_octave
%|function out = ir_is_octave
%|
%| determine if this is octave running
%| inspired by http://wiki.octave.org/Compatibility
%|
%| out is true or false
%|
%| Jeff Fessler

persistent x
if (isempty(x))
	x = exist('OCTAVE_VERSION', 'builtin');
	x = (x == 5); % 5 signifies a built in function
end
	out = x;
end

%tmp = version;
%out = streq(tmp, '3.6.3');
