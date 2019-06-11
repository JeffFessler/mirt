 function yn = has_mex_jf
%function yn = has_mex_jf
% determine whether this matlab session has access to .mex files needed
% for advanced features of this toolbox.

yn = which('dtft_mex'); % seems to force a path check ??

yn =	...
	exist('jf_mex') == 3 && ...
	exist('wtfmex') == 3 && ...
	exist('f3d_mex') == 3 && ...
	exist('dtft_mex') == 3 && ...
	exist('rotmex') == 3;

persistent Mexwarned
if ~isvar('Mexwarned') || isempty(Mexwarned)
	Mexwarned = false;
end

if ~yn && ~Mexwarned
	if ispc
		warn('your operating system came from a monopoly.')
		warn('mex files are not available for this system.')
		warn('upgrade to unix to get the full features!')
	else
		warn('mex files not found: are your paths correct?')
	end
	Mexwarned = true;
end
