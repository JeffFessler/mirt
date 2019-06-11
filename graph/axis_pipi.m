 function axis_pipi(varargin)
%function axis_pipi(['set',] xo, yo, omx, omy)
%|
%| label x and y axes from -pi to pi for 2D DSFT
%| (if a 3D plot, then also label z axis from -pi to pi for 3D DSFT)

if nargin <= 1
	is3d = numel(axis) == 6;
	% trick: axes labels must *precede* the symbol font below
	xlabelf '$\omx$'
	ylabelf '$\omy$'
	if is3d
		zlabelf '$\omz$'
	end
	xaxis_pi '-p 0 p'
	yaxis_pi '-p 0 p'
	if is3d
		zaxis_pi '-p 0 p'
	end

return
	if nargin == 1
		if streq(varargin{1}, 'set')
			axis([-pi pi -pi pi])
		else
			error 'unknown option'
		end
	end
%	axis square
	xtick([-pi 0 pi])
	ytick([-pi 0 pi])
	set(gca, 'xticklabel', '-p | 0 | p', 'fontname', 'symbol')
	set(gca, 'yticklabel', '-p | 0 | p', 'fontname', 'symbol')
return
end

warning 'the old style is obsolete!'

% fix: use varargin hereafter...
if ~isvar('xo') || isempty(xo), xo = 0.05; end
if ~isvar('yo') || isempty(yo), yo = 0.05; end
if ~isvar('omx'), omx = '\omega_X'; end
if ~isvar('omy'), omy = '\omega_Y'; end

t = axis;
if (strcmp(get(gca, 'XScale'), 'log'))
	t = log10(t);
	xo = 10^(t(1) - (t(2)-t(1))*xo);
else
	xo = t(1) - (t(2)-t(1))*xo;
end

t = axis;
if (strcmp(get(gca, 'YScale'), 'log'))
	t = log10(t);
	yo = 10^(t(3) - (t(4)-t(3))*yo);
else
	yo = t(3) - (t(4)-t(3))*yo;
end

axis(pi * [-1 1 -1 1])
set(gca, 'xtick', [-pi -pi/2 0 pi/2 pi])
set(gca, 'xticklabel', [])
%set(gca, 'xticklabel', ['-pi  '; '-pi/2'; '  0  '; 'pi/2 '; 'pi   '])
set(text(pi, yo, '\pi'), 'horizontalalignment', 'center')
set(text(-pi, yo, '-\pi'), 'horizontalalignment', 'center')
set(text(pi/2, yo, '\pi/2'), 'horizontalalignment', 'center')
set(text(-pi/2, yo, '-\pi/2'), 'horizontalalignment', 'center')
%set(text(0, yo, '0'), 'horizontalalign', 'center')

if ~isempty(omx)
	set(text(0, yo, omx), 'horizontalalignment', 'center')
end

set(gca, 'ytick', [-pi -pi/2 0 pi/2 pi])
set(gca, 'yticklabel', [])
rtext(xo, pi, '\pi')
rtext(xo, -pi, '-\pi')
rtext(xo, pi/2, '\pi/2')
rtext(xo, -pi/2, '-\pi/2')

%rtext(xo, 0, '0')
if ~isempty(omy)
	rtext(xo, 0, omy)
end

function rtext(x, y, s)
	set(text(x, y, s), ...
	'verticalalignment', 'middle', ...
	'horizontalalignment', 'right');
