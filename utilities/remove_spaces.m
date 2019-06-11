  function arg = remove_spaces(arg)
%|function arg = remove_spaces(arg)
%| replace extra spaces at ends of matrix string array with zeros.
%| also add extra '0' at end to be sure to null terminate string
%| Copyright May 2000, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

if streq(arg, 'test'), remove_spaces_test, clear arg, return, end

for ii=1:size(arg,1)
	jj = size(arg,2);
	while (arg(ii,jj) == ' ')
		arg(ii,jj) = 0;
		jj = jj - 1;
		if jj < 1
			warn('bug in %s', mfilename)
			keyboard
		end
	end
end
arg(:,end+1) = 0;

arg(arg(:,1) == 10, :) = []; % remove blank lines too

function remove_spaces_test
str = strvcat('a', 'bb ', 'c c ');
out = remove_spaces(str);
tmp = char([97 0 0 0 0; 98 98 0 0 0; 99 32 99 0 0]);
jf_equal(out, tmp)
