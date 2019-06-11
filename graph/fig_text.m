  function hh = fig_text(varargin)
%|function hh = fig_text(...)
%| add text to figure, using coordinates relative to entire window
%| options:
%|	'-date'		prepend date
%|	'-tex'		tex parse string
%|	x,y		coordinates relative to [0 0 1 1]
%|	{args}		style arguments for text() command
%| jeff fessler

comment = '';
thedate = '';
interpret = 'none';
textarg = {'fontsize', 11};
x = 0.01;
y = 0.01;
while length(varargin)
	arg = varargin{1};
	if ischar(arg)
		if streq(arg, '-date')
			thedate = [mydate ' '];
		elseif streq(arg, '-tex')
			interpret = 'tex';
		else
			comment = arg;
		end
		varargin = {varargin{2:end}};

	elseif isnumeric(arg)
		if length(varargin) < 2, error 'need x, y', end
		x = varargin{1};
		y = varargin{2};
		varargin = {varargin{3:end}};

	elseif iscell(arg)
		textarg = {textarg{:}, arg{:}};
		varargin = {varargin{2:end}};

	else
		help(mfilename), error 'bad args'
	end
end

hh = gca;

hold on
h = axes('position', [0 0 1 1]); axis off
h = text(x, y, [thedate comment], 'interpreter', interpret, textarg{:});
hold off

if nargout
	hh = h;
else
	axes(hh) % reset axes to the one before the text added
	clear hh
end

function s = mydate
%s = datestr(now, 'yyyy-mm-dd HH:MM:SS');
%s = datestr(now, 'yy-mm-dd HH:MM:SS');
t1 = datestr(now, 'yy');
t2 = datestr(now, 'mm');
t3 = datestr(now, 'dd');
t4 = datestr(now, 'HH:MM:SS');
s = sprintf('%s-%s-%s %s', t1, t2, t3, t4);
