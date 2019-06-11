  function out = jf_show_iter(varargin)
%|function out = jf_show_iter([options])
%|
%| function handle to serve as 'userfun' in iterative algorithms
%| to show current iterate
%|
%| call it first to initialize with options:
%|	'mask'	support mask for embed()
%|	'name'	variable name.  default: 'x'
%|	'clim'	color limits for im().  default: []
%|	'draw'	call drawnow()?  default: true
%|	'pause'	call pause() after each display?
%|
%| Copyright 2010-03-10, Jeff Fessler, University of Michigan

persistent arg
if ~isvar('arg') || isempty(arg)
	arg.mask = [];
	arg.name = 'x';
	arg.clim = [];
	arg.draw = true;
	arg.pause = false;
end

% trick: called with (x, ...)
if numel(varargin) >= 1 && isnumeric(varargin{1})
	x = varargin{1};
%	pr size(x)
	varargin = {varargin{2:end}};
end

% trick: called with (x, iter, ...)
if numel(varargin) >= 1 && isnumeric(varargin{1})
	iter = varargin{1};
	pr iter
	varargin = {varargin{2:end}};
end

if numel(varargin) >= 1
	if ischar(varargin{1})
		printm 'setup'
		arg = vararg_pair(arg, varargin);
		out = 0;
	else
		fail 'bad usage'
	end
return
end

if ~isvar('x')
	x = evalin('caller', arg.name);
end
if isempty(arg.clim)
	im(embed(x, arg.mask))
else
	im(embed(x, arg.mask), arg.clim)
end

if arg.draw
	drawnow
end
if arg.pause
	pause
end
out = 0;
