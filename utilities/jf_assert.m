 function jf_assert(varargin)
%function jf_assert(command)
% verify that the command (evaluated within caller) returns true.
% if not, print error message.

if nargin < 1, ir_usage, end
if nargin == 1 && streq(varargin{1}, 'test'), jf_assert_test, return, end

arg = [varargin{:}]; % handle cases with spaces like 'jf_assert x == y'

[name line] = caller_name;
if isempty(name)
	str = '';
else
	str = sprintf(' at %d in "%s"', line, name);
end

try
	tmp = evalin('caller', arg);
catch
	error(['%s was unable to evaluate "%s"' str], mfilename, arg)
end

% note: use isequal, not issame!
if isscalar(tmp) && islogical(tmp)
	if ~tmp
		fail(['jf_assert of "%s" was untrue' str], arg)
%		dbup
%		dbstack
%		keyboard
	end
else
tmp
	whos
	error(['jf_assert of "%s" did not return logical scalar' str], arg)
end

function jf_assert_test
jf_assert 7 == 7
