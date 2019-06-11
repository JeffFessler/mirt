  function ob = Gbspline1(varargin)
%|function ob = Gbspline1(options)
%|
%| Construct Gbspline1 object for 1D B-spline interpolation
%|
%| See Gbspline1_test.m for example usage.
%|
%| in
%|
%| options
%|	'type'	char		'val2coef' convert signal values to coefficients
%|				'
%|				'coef2val' | 'synth' conv
%|	'mask'	size(image)	logical array of object support.
%|	'chat'		verbose printing of debug messages
%|	'psf'		point spread function (aka impulse response)
%|	'type'		type of bspline1:
%|				'conv,same'	usual case (default)
%|				'fft,same'	(periodic end conditions)
%|				todo: allow replicated end conditions!
%|				'imfilter,same'	(requires image toolbox)
%|				'imfilter,circ'	(requires image toolbox)
%|	'imfilter_options'	options to 'imfilter' (if used)
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|			nd = prod(size(mask)) for 'conv,same' type
%|
%| Copyright 2011-11-13, Jeff Fessler, University of Michigan

fail 'not done'

if nargin == 1 && streq(mask, 'test'), Gbspline1_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

% option defaults
arg.type = 'todo';
arg.mask = mask;
arg.chat = 0;

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

does_many = false;
switch arg.type
case 'todo'
	arg.odim =
	handle_forw = @(arg,x)
	handle_back = @(arg,y)
	does_many = true;

otherwise
	error 'unknown bspline1 type'
end

ob = fatrix2('idim', arg.odim, 'arg', arg, 'does_many', does_many, ...
	'forw', handle_forw, 'back', handle_back);
