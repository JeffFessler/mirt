  function hs = jf_add_slider(varargin)
%|function hs = jf_add_slider([options])
%|
%| add a slider control to the current figure to enable interaction
%|
%| option
%|	'callback'	@ (value,data{:})	function to call with
%|						slider value [0,1] and data
%|	'pos'		[l b w h]	slider position
%|	'data'		{}		data to be passed to callback
%|	'value'		0 <= x <= 1	initial slider value, default 0
%|
%| Copyright 2010-1-11, Jeff Fessler, University of Michigan

if ~nargin && ~nargout, help(mfilename), error(mfilename), end
if nargin && streq(varargin{1}, 'test'), jf_add_slider_test, return, end

arg.callback = @jf_add_slider_call;
arg.pos = [0.1, 0.0, 0.8, 0.03];
arg.data = {};
arg.name = '';
arg.value = 0;
arg.min = [];
arg.max = [];
arg.style = 'slider';
arg.sliderstep = [0.01 0.10];
arg.string = 'button';
arg = vararg_pair(arg, varargin);

args = {};
if ~isempty(arg.min)
	args = {args{:}, 'min', arg.min};
end
if ~isempty(arg.max)
	args = {args{:}, 'max', arg.max};
end

if isempty(arg.callback)
	arg.callback = @jf_add_slider_call_null;
end

% construct gui components
if ir_is_octave % && exist('uicontrol') ~= 2
	warn 'no uicontrol in octave'
	% todo: draw box or something?
else
	hs = uicontrol('style', arg.style, 'string', arg.string, ...
		'units', 'normalized', 'position', arg.pos, ...
		'value', arg.value, ...
		args{:}, ...
		'sliderstep', arg.sliderstep, ...
		'callback', {@jf_add_slider_callback});
end

if ~isempty(arg.name)
	set(gcf, 'name', arg.name) % assign GUI a name for window title
end

if ~nargout
	clear hs
end


% callbacks: these automatically have access to component handles
% and initialized data because they are nested at a lower level.

% slider callback
function jf_add_slider_callback(source, eventdata) 
	val = get(source, 'value');
	arg.callback(val, arg.data{:});
end % jf_add_slider_callback

end % jf_add_slider


function jf_add_slider_call(value)
pr value
end % jf_add_slider_call

function jf_add_slider_call_null(value)
end % jf_add_slider_call_null


function jf_add_slider_test
jf_add_slider_test_call(0.5);
if im
	jf_add_slider('callback', @jf_add_slider_test_call, 'data', {gca}, ...
		'value', 0.5)
%	jf_add_slider_test_call(0);
	jf_add_slider('callback', @jf_add_slider_test_call, 'data', {gca}, ...
		'pos', [0.1 0.04 0.8 0.03], ...
		'value', 0, 'style', 'togglebutton', 'string', 'push me')
end
end % jf_add_slider_test


function jf_add_slider_test_call(value, h)
if isvar('h')
	axes(h)
end
if im
	plot([0 value], [0 value], '-')
	axis([0 1 0 1])
end
end % jf_add_slider_test_call
