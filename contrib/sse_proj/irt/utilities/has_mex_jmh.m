function yn = has_mex_jmh
% for advanced features of this toolbox.

yn = which('jmh_sse_mex'); % seems to force a path check ??

yn =	...
	exist('jmh_sse_mex') == 3;

persistent Mexwarned
if isempty(Mexwarned)
	Mexwarned = 0;
end

if ~yn & ~Mexwarned
  warning('mex files not found: are your paths correct?')
  Mexwarned = 1;
end
