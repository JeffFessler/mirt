 function yn = jf_has_hct2(pn)
%function yn = jf_has_hct2(pn)
% determine if this matlab session has access to the binary "hct2" executable,
% for internal UM use only.

yn = 0;

if ~isunix, return, end

try
	% unix returns nonzero on failure
	if unix('which hct2 >& /dev/null')
		return
	end

	yn = 1;
	return

catch
	printm('the last error was "%s"', lasterr)
	return
end

error 'bug'
