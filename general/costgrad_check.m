  function costgrad_check(x, data, costgrad, varargin)
%|function costgrad_check(x, data, costgrad, varargin)
%|
%| Check for consistency between a cost function and its gradient.
%|
%| in
%|	x	[np 1]		point at which to evaluate cost and gradient
%|	data	{cell}		whatever data is needed for the cost function
%|	costgrad @()		function returning cost function and gradient:
%|					[cost grad] = costgrad(x, data)
%|
%| option
%|	'step'			user-selected initial step (default 0.01)
%|
%| 2010-10-28, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test')
	run_mfile_local('costgrad_check_test'), return
end
if nargin < 3, help(mfilename), error(mfilename), end
if ~isa(costgrad, 'function_handle'), error 'costgrad not function handle?', end

arg.step = 0.01; % dumb user default

arg = vararg_pair(arg, varargin);
arg.isave = iter_saver(arg.isave, arg.niter);

% gradient of cost function
[cost0 grad] = costgrad(x, data);

dir = zeros(size(x));
jj = 1;
dir(jj) = arg.step;
x1 = x + dir;

cost1 = costgrad(x1, data);

der = (cost1 - cost0) / arg.step;
pr [der grad(jj)]


function costgrad_check_text
fail todo
