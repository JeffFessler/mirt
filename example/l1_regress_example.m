%% l1_regress_example.m
%
% Example of $l_1$ regression for robust estimation
% $\min_x \| y - A x \|_1$
% except that $l_1$ is approximated by a hyperbola (corner rounding).
%
% 2005-4-24, Jeff Fessler, University of Michigan

%% Example data with one outlier point

ti = [0:10]'; % sample locations
rng(0)
yi = 3 + 2 * ti;
yi = yi + 0.2 * randn(size(yi));
yi(9) = yi(1); % outlier
if im
	plot(ti, yi, '-o')
end

%% Iterative reweighted LS algorithm (Huber's method)

if ~isvar('x1') || 0, printm 'l1 regression'
	f.niter = 10;
	f.delta = 0.01;
	% must choose 'delta' small enough to reject outliers, but not
	% so small that convergence is too slow.  play around with it!
	[x1s, x2] = l1_regress_fun(ti, yi, ... % data
		'niter', f.niter, 'delta', f.delta, ... % parameters
		'linear', false); % use affine
	x1 = x1s(:,end);
end

%% ADMM version

if ~isvar('x1_admm') || 1, printm 'l1 regression via admm'
	A = [ones(size(ti)) ti];
	x1s_admm = l1_regress_admm1(yi, A, 'niter', f.niter, ...
		'x0', x1, ...
		'x0', [0; 0], ...
		'x0', [], ...
		'rho', 2^-3, ... % trick: had to tune :(
		'isave', 'all');
	x1_admm = x1s_admm(:,end);
end

%% Plot parameters vs iteration and robust regression lines
% The parameters appear to converge (essentially) within 8-10 iterations.

if im
	clf
	subplot(211)
	plot(0:f.niter, x1s', '-o', 0:f.niter, x1s_admm', '-+')
	xlabel 'iteration'
	axisy(0,4)
	ylabel 'estimate'

	tt = linspace(min(ti), max(ti), 101);

	subplot(212)
	plot(ti, yi, 'bo', ...
		tt, x1(1) + x1(2) * tt, 'g-', ...
		tt, x1_admm(1) + x1_admm(2) * tt, 'm:', ...
		tt, x2(1) + x2(2) * tt, 'r--')
	ir_legend({'data', '$l_1$ regression via Huber', ...
		'$l_1$ via admm', '$l_2$ regression'})
	xlabelf('$t_i$')
	ylabelf('$y_i$')
end
