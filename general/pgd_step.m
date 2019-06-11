  function [xs, info] = pgd_step(x, data, costgrad, varargin)
%|function [xs, info] = pgd_step(x, data, costgrad, varargin)
%|
%| Generic unconstrained minimization
%| via preconditioned gradient descent algorithm with diminishing step size,
%| e.g., reducing step size by half until cost function decreases.
%| (This method is practical only when cost function is computed fairly fast.)
%|
%| in
%|	x	[np 1]		initial estimate
%|	data	{cell}		whatever data is needed for the cost function
%|	costgrad @()		function returning cost function and gradient:
%|					[cost grad] = costgrad(x, data)
%|
%| option
%|	'precon' [np np]	preconditioner (default 1)
%|	'niter'			# total iterations (default 1)
%|	'step'			user-selected initial step (default 1)
%|	'step_method'		'der2' from 2nd derivative (default) or 'user'
%|	'step_fraction'		decrease step by this each try (default 0.5)
%|	'max_num_step'		max # of tries at decreasing step (default 99)
%|	'chat'			verbosity (default 0)
%|
%| out
%|	xs	[np niter]	estimates each iteration
%|	info	[niter 3]	step, time, cost
%|
%| Copyright 2005-1-10, Anastasia Yendiki, University of Michigan
%| 2010-10 modified by Jeff Fessler and Amanda Funai

if nargin == 1 && streq(x, 'test')
	run_mfile_local('pgd_step_test'), return
end
if nargin < 3, help(mfilename), error(mfilename), end
if ~isa(costgrad, 'function_handle'), error 'costgrad not function handle?', end

arg.precon = 1; % default precondioner
arg.niter = 1;
arg.step_method = 'der2'; % default is to try 2nd derivative
arg.step = 1; % dumb user default
arg.step_fraction = 0.5;
arg.chat = 0;
arg.eps_fraction = 0.01;
arg.eps = '';
arg.eps_lower = 0.001; % see below
arg.eps_upper = 0.1;
arg.eps_warn = false; % not sure that eps_upper _lower are useful
arg.max_num_step = 99; % at most this many halving attempts
arg.isave = [];

arg = vararg_pair(arg, varargin);
arg.isave = iter_saver(arg.isave, arg.niter);

xs = zeros(numel(x), length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

info = zeros(arg.niter, 3);

cpu etic
for iter=1:arg.niter
	ticker(mfilename, iter, arg.niter)

	% gradient of cost function
	[cost0 grad] = costgrad(x, data);

	% preconditioned gradient
	if isa(arg.precon, 'function_handle')
		pregrad = arg.precon(x) * grad;
	else
		pregrad = arg.precon * grad;
	end

	% search direction
	ddir = -reshape(pregrad, size(x));

	if all(ddir == 0)
		warn('done at iteration %d', iter)
		xs(:,arg.isave >= iter) = repmat(x, [1 sum(arg.isave >= iter)]);
		return
	end

	switch arg.step_method
	case 'der2' % determine step size by numerical 2nd derivative
		if isempty(arg.eps)
			eps0 = max(abs(x)) / max(abs(ddir)) * arg.eps_fraction;
		else
			eps0 = arg.eps;
		end
		% 1st and 2nd derivatives of f(step) = cost(x + step * ddir)
		[cost_eps grad_eps] = costgrad(x + eps0 * ddir, data);
		der1 = real(sum(col(conj(ddir) .* grad)));
		der2 = (real(sum(col(conj(ddir) .* grad_eps))) - der1) / eps0;
		if der1 < 0 && der2 > 0
			step1 = -der1 / der2; % newton around step=0
			if arg.eps_warn && eps0 / step1 < arg.eps_lower
				warn('eps=%g small vs step=%g', eps0, step1)
			end
			if arg.eps_warn && eps0 / step1 > arg.eps_upper
				warn('eps=%g large vs step=%g', eps0, step1)
			end
		else
			warn('der1=%g, der2=%g', der1, der2)
			warn '2nd derivative step failed; revert to user value'
			step1 = arg.step; % revert to user-selected value
		end
	case 'user'
		step1 = arg.step;
	otherwise
		fail('unknown step_method "%s"', arg.step_method)
	end

	steps = zeros(arg.max_num_step, 1);
	costs = zeros(arg.max_num_step, 1);
	for is=1:arg.max_num_step
		step = step1 * (arg.step_fraction)^(is-1);
		xnew = x + step * ddir;
		costn = costgrad(xnew, data);
		steps(is) = step;
		costs(is) = costn;
		if (costn <= cost0) % done if cost decreased
			is = 0;
			break;
		end
		if arg.chat && is > 1
			printm 'warn: decreasing step again'
		end
	end

	% if we could not descend, show diagnostic plot
	if is ~= 0
		plot(steps, costs, '.-', 0, cost0, '*')
		xlabel 'step', ylabel 'cost'
		fail('error: too many half steps')
	end

	% update
	x = xnew;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end

	info(iter,:) = [step cpu('etoc') costn];
end
