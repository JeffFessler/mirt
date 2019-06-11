 function ir_stem_spectrum(x, y, varargin)
%function ir_stem_spectrum(x, y, varargin)
%| 
%| default spectrum style for engn 100

arg.thresh = 0;
arg.color = [0 0.5 0];
arg = vararg_pair(arg, varargin);

k = y > max(y(:)) * arg.thresh;

stem(x(k), y(k), ...
	'color', arg.color, ...
	'linewidth', 3.0, ...
	'marker', 'none')
ylabelf('amplitude')
