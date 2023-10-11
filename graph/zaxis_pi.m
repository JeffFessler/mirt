 function zaxis_pi(varargin)
%function zaxis_pi(varargin)
%| label z axis with various forms of "pi"
%| the argument can be a string with p's in it, or fractions of pi:
%| [0 1/2 1] or '0 p/2 p' -> [0 pi/2 pi]
%| [-1 0 1] or '-p 0 p' -> [-pi 0 pi]
%| etc.
%|
%| There is one catch: this changes the axes font, so subsequent
%| calls to title or xlabel or ylabel will have the wrong font.
%| so title, xlabel, ylabel should be done *before* calling this routine.
%|
%| Jeff Fessler

arg.do_tex = 1;

if length(varargin) == 0
	ticks = '0 p';
elseif length(varargin) == 1
	ticks = varargin{1};
else
	fail 'only one arg allowed'
end

if ischar(ticks)
	str = ticks;
	str = strrep(str, ' ', ' | ');
	str = strrep(str, '*', ''); % we don't need the "*" in label
	ticks = strrep(ticks, '2p', '2*p');
	ticks = strrep(ticks, '3p', '3*p');
	ticks = strrep(ticks, '4p', '4*p');
	ticks = strrep(ticks, '5p', '5*p');
	ticks = strrep(ticks, 'p', 'pi');
	ticks = eval(['[' ticks ']']);

else

	if same(ticks, [0])
		str = '0';
	elseif same(ticks, [0 1])
		str = '0 | p';
	elseif same(ticks, [0 1/2 1])
		str = '0 | p/2 | p';
	elseif same(ticks, [-1 0 1])
		str = '-p | 0 | p';
	elseif same(ticks, [0 1 2])
		str = '0 | p | 2p';
	else
		fail 'this ticks not done'
	end

end

% here is the main part
zlim([min(ticks), max(ticks)])
ztick(ticks)

%arg.do_tex = streq(version('-release'), '2014b'); % supports TeX in ticks

if ~ir_is_octave && arg.do_tex
	str = strrep(str, 'p', '\pi');
	str = strrep(str, ' | ', ' ');
	if str(end) ~= ' ', str(end+1) = ' '; end % space at end
	tmp = strfind(str, ' '); % find all spaces
	ntick = numel(tmp);
	stick = cell(ntick,1);
	tmp = [0 tmp];
	for ii=1:ntick
		stick{ii} = str((tmp(ii)+1):(tmp(ii+1)-1));
	end
	h = gca;
	h.ZTickLabel = stick;
return
end

set(gca, 'zticklabel', str)
ir_set_gca_fontname('symbol')

if ir_is_octave
	persistent warned
	if isempty(warned)
		warned = true;
		warn 'may not print properly'
	end
end


function is = same(x,y)
if length(x) ~= length(y)
	is = 0;
	return
end
is = all(x == y);
