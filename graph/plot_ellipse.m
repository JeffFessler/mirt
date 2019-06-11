  function [xo, yo] = plot_ellipse(cx, cy, rx, ry, theta, varargin)
%|function [xo, yo] = plot_ellipse(cx, cy, rx, ry, theta, [options])
%| plot an ellipse
%| option
%|	'n'		# of points for half of ellipse.  default: 301
%|	'c'		line type.  default: 'b-'
%|	'hold'	1|0	if 1, then add to current plot.  default: 0

if nargin == 1 && streq(cx, 'test'), plot_ellipse_test, return, end
if nargin < 5, help(mfilename), error(mfilename), end

arg.n = 301;
arg.c = 'b-';
arg.hold = false;
arg = vararg_pair(arg, varargin);

xe = linspace(-rx, rx, arg.n)';
yp = ry * sqrt(1 - (xe/rx).^2);
ym = -yp;
xp = cx + (cos(theta) * xe - sin(theta) * yp);
yp = cy + (sin(theta) * xe + cos(theta) * yp);
xm = cx + (cos(theta) * xe - sin(theta) * ym);
ym = cy + (sin(theta) * xe + cos(theta) * ym);

xo = [xp; flipud(xm)];
yo = [yp; flipud(ym)];

if ~nargout
	if im
		if arg.hold, hold on, end
		plot(xo, yo, arg.c)
		if arg.hold, hold off, end
	end
	clear xo yo
end

function plot_ellipse_test
plot_ellipse(5, 10, 10, 20, 1)
[xo yo] = plot_ellipse(5, 10, 10, 20, 1);
