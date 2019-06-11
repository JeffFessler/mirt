 function jpg_write(x, file, varargin)
%function jpg_write(x, file, varargin)

arg.clim = [];
arg.qual = 99;

arg = vararg_pair(arg, varargin);

if isempty(arg.clim)
	arg.clim = minmax(x);
end

x = floor(255*(x-arg.clim(1))/(arg.clim(2)-arg.clim(1)));
x = max(x,0);
x = min(x,255);
x = uint8(x);
x = x.';

file = file_check(file, '.jpg');

if ~isempty(file)
	imwrite(x, file, 'quality', arg.qual)
end


function file = file_check(base, suff)

keepask = 1;
while (keepask)
	keepask = 0;
	file = [base suff];

	if exist(file, 'file')
		prompt = sprintf('overwrite figure "%s"? y/a/[n]: ', file);
	else
		prompt = sprintf('print new figure "%s"? y/n: ', file);
	end

	ans = input(prompt, 's');
	% alternate filename
	if streq(ans, 'a')
		t = sprintf('enter alternate basename (no %s): ', suff);
		base = input(t, 's');
		keepask = 1;
	end
end

if ~isempty(ans) && ans(1) == 'y'
	return
end
file = [];
