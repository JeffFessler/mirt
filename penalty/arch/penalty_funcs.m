 function [potk, wpot, dpot] = penalty_funcs(type, beta, delta)
%function [potk, wpot, dpot] = penalty_funcs(type, beta, delta)
%
% Define roughness penalty functions via inline functions.
%
% If beta is numeric, then the penalty will have the form
%	R(x) = sum_k beta * potential([Cx]_k, delta)
% If beta is empty, then the penalty will have the form
%	R(x) = sum_k wt_k * potential([Cx]_k, delta)
% in:
%	type		quad or huber or lange3 or ...
%	beta and delta should be numeric values (full precision!)
% out:
%	potk(wt, C*x)	potential function (inline)
%	wpot(wt, C*x)	potential 'weights' (ala half-quad. curvatures)
%	dpot(wt, C*x)	potential derivative
%
% Copyright 2000-5, Jeff Fessler, The University of Michigan
% 2002-2-20 changed so that funcitons have TWO arguments (allow weighting)
% 2005-2-14 discovered that 'hyper' was missing sqrt()

if nargin < 1, help(mfilename), error(mfilename), end

% default is to plot, thereby testing
if streq(type, 'test'), penalty_funcs_test, return, end

if ~isvar('beta'), beta = []; end


%
% quadratic potential function
%
if streq(type, 'quad')
	if isempty(beta)
		potk = '(abs(t).^2)/2';
		wpot = 'ones(size(t))';
		dpot = 't';
	else
		warning 'this usage is obsolete'
		potk = sprintf('%19.18e .* (abs(t).^2)/2', beta);
		wpot = sprintf('%19.18e .* ones(size(t))/2', beta);
		dpot = sprintf('%19.18e .* t', beta);
	end


%
% huber potential function
%
elseif streq(type, 'huber')
	if isempty(beta)
		potk = sprintf('huber_pot(t,%19.18e)', delta);
		wpot = sprintf('huber_wpot(t,%19.18e)', delta);
		dpot = sprintf('huber_dpot(t,%19.18e)', delta);
	else
		warning 'this usage is obsolete'
		potk = sprintf('%19.18e .* huber_pot(t,%19.18e)', beta, delta);
		wpot = sprintf('%19.18e .* huber_wpot(t,%19.18e)', beta, delta);
		dpot = sprintf('%19.18e .* huber_dpot(t,%19.18e)', beta, delta);
	end


%
% this is a non-convex penalty that was my initial incorrect
% attempt at a hyperbola penalty
% but actually it was the cauchy prior: d^2 / 2 * log(1 + (t/d)^2)
%
elseif streq(type, 'cauchy')
	if isempty(beta)
%		potk = '%19.18e .* (sqrt(1 + abs(t/%19.18e).^2) - 1)';
		potk = '%19.18e / 2 .* log(1 + abs(t/%19.18e).^2)';
		potk = sprintf(potk, delta.^2, delta);
		wpot = sprintf('1 ./ (1 + abs(t/%19.18e).^2)', delta);
		dpot = sprintf('t ./ (1 + abs(t/%19.18e).^2)', delta);
	else
		warning 'this usage is obsolete'
		potk = '%19.18e .* (sqrt(1 + abs(t/%19.18e).^2) - 1)';
		potk = sprintf(potk, beta*delta^2, delta);
		wpot = sprintf('%19.18e ./ (1 + abs(t/%19.18e).^2)', beta, delta);
		dpot = '%19.18e .* t ./ (1 + abs(t/%19.18e).^2)';
		dpot = sprintf(dpot, beta, delta);
	end


elseif streq(type, 'hyper2')
	if isempty(beta)
		potk = '%19.18e .* (sqrt(1 + abs(t/%19.18e).^2) - 1)';
		potk = sprintf(potk, delta.^2, delta);
		wpot = sprintf('1 ./ sqrt(1 + abs(t/%19.18e).^2)', delta);
		dpot = sprintf('t ./ sqrt(1 + abs(t/%19.18e).^2)', delta);
	else
		warning 'this usage is obsolete'
		potk = '%19.18e .* (sqrt(1 + abs(t/%19.18e).^2) - 1)';
		potk = sprintf(potk, beta*delta^2, delta);
		wpot = sprintf('%19.18e ./ sqrt(1 + abs(t/%19.18e).^2)', beta, delta);
		dpot = '%19.18e .* t ./ sqrt(1 + abs(t/%19.18e).^2)';
		dpot = sprintf(dpot, beta, delta);
	end

elseif streq(type, 'hyper')
	'the original "hyper" potential was wrong.'
	'choose "cauchy" to get the original (non-convex!) hyper effect.'
	'choose "hyper2" for the corrected version'
	error ''


%
% Lange3 penalty
%
elseif streq(type, 'lange3')
	if isempty(beta)
		potk = '%19.18e .* (abs(t)/%19.18e - log(1+abs(t)/%19.18e))';
		potk = sprintf(potk, delta^2, delta, delta);
		wpot = sprintf('1 ./ (1 + abs(t)/%19.18e)', delta);
		dpot = sprintf('t ./ (1 + abs(t)/%19.18e)', delta);
	else
		warning 'this usage is obsolete'
		potk = '%19.18e .* (abs(t)/%19.18e - log(1+abs(t)/%19.18e))';
		potk = sprintf(potk, beta*delta^2, delta, delta);
		wpot = sprintf('%19.18e ./ (1 + abs(t)/%19.18e)', beta, delta);
		dpot = '%19.18e .* t ./ (1 + abs(t)/%19.18e)';
		dpot = sprintf(dpot, beta, delta);
	end


else
	error(sprintf('Unknown type "%s"', type))
end

potk = inline(['wt .* ' potk], 'wt', 't');
wpot = inline(['wt .* ' wpot], 'wt', 't');
dpot = inline(['wt .* ' dpot], 'wt', 't');


function penalty_funcs_test
plist = {'quad', 'huber', 'lange3', 'cauchy', 'hyper2'};
t = linspace(-20,20,101)';
del = 5;
for ii=1:length(plist)
	type = plist{ii};
	[potk wpot dpot] = penalty_funcs(type, [], del);
	pp(:,ii) = potk(1, t);
	pw(:,ii) = wpot(1, t);
	pd(:,ii) = dpot(1, t);
end
if im
	clf
	subplot(311)
	plot(t, pp), title potk, legend(plist)
	subplot(312)
	plot(t, pw), title wpot, legend(plist)
	subplot(313)
	plot(t, pd), title dpot, legend(plist)
	axisy(-5.5,5.5)
end
