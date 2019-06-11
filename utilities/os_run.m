  function out = os_run(str)
%|function out = os_run(str)
%| call OS (unix only of course), check for error, optionally return output

if nargin < 1, ir_usage, end
if streq(str, 'test'), os_run_test, return, end

[s out1] = unix(str);
if s
	fail('unix call failed:\n%s', str)
end

if nargout
	out = out1;
end

function os_run_test
printm 'os_run test'
if ~isunix
	warn 'os_run works only on unix'
	return
end

out = os_run('echo 1+2 | bc');
jf_equal(out, sprintf('3\n'))
