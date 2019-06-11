 function yn = has_aspire
%function yn = has_aspire
% determine whether this matlab session has access to the Aspire
% iterative reconstruction executables.
% those executables are not needed for most purposes, except for my
% own testing of consistency between Matlab and C implementations.

yn = 0;

if ~isunix, return, end

try
	% unix returns nonzero on failure
	if unix('which op >& /dev/null') ...
		|| unix('which wt >& /dev/null') ...
		|| unix('which i >& /dev/null')
		return
	end

	yn = 1;
	return

catch
	printf('the last error was "%s"', lasterr)
	return
end

error 'bug'
