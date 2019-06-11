% ir_fista_test
% test FISTA using the quadratic programming example from O'Donoghue 2015
% min_x 1/2 x' Q x + q' x sub to a <= x <= b
% 2017 Jeff Fessler with corrections by Donghwan Kim

if ~isvar('step')
	rng(0)
	np = 500;
	Q = randn(np,np);
	Q = Q' * Q;
%	pr cond(Q)
	tmp_eps = 10^-4;
%	tmp_eps = 10^-7; % constrained
	Q = (1-tmp_eps) * Q/max(eig(Q)) + tmp_eps * eye(np,np);
%	pr cond(Q) % about 1/eps

	a = -ones(np,1); % lower bound
	b = ones(np,1); % upper bound
	q = 0.005 * randn(np,1); % todo: why
	step = 1 / eigs(Q,1)

% todo: check broadcast
	f.cost_jf = @(x) 0.5 * sum(x .* (Q * x), 1, 'double') - q' * x;
	f.cost = @(x) 0.5 * sum(x .* (Q * x - repmat(q,1,size(x,2))), 1, 'double'); % DK
	f.grad = @(x) Q * x - q;

	if 0 % constrained
		f.prox = @(x) max(min(x,b), a);
	else % unconstrained
		f.prox = @(x) x;
		xhat = Q \ q;
		f.min = f.cost(xhat);
	end

	x0 = zeros(np,1); % initial guess
%	cost0 = f.cost(x0);
	f.niter = 2000;
end

%% GP
if ~isvar('xgp'), printm 'gradient projection'
	xgp = ir_fista(x0, 'gradfun', f.grad, 'step', step, 'restart', 0, ...
		'proxfun', f.prox, ...
		'niter', f.niter, 'isave', 'all', 'momentum', 0);
	cost.gp = f.cost(xgp);
end

%% FPGM
if ~isvar('xfg'), printm 'fista'
	xfg = ir_fista(x0, 'gradfun', f.grad, 'step', step, 'restart', 0, ...
		'proxfun', f.prox, ...
		'niter', f.niter, 'isave', 'all');
	cost.fg = f.cost(xfg);
end

%% FPGM + restart
% todo: why is this restart not helping?
if ~isvar('xar'), printm 'adaptive restart'
	xar = ir_fista(x0, 'gradfun', f.grad, 'step', step, 'restart', 1, ...
		'proxfun', f.prox, ...
		'niter', f.niter, 'isave', 'all', 'chat', 1);
	cost.ar = f.cost(xar);
end

if ~isvar('xar90'), printm 'adaptive restart with 90'
	xar90 = ir_fista(x0, 'gradfun', f.grad, 'step', step, 'restart', 1, ...
		'proxfun', f.prox, ...
		'alpha_restart', 0, ...
		'niter', f.niter, 'isave', 'all', 'chat', 1);
	cost.ar90 = f.cost(xar90);
end

%% todo: this is viable only if ir_fista is working well!
if ~isvar('xinf'), printm 'xinf'
	xinf = ir_fista(xar(:,end), ...
		'gradfun', f.grad, 'step', step, 'restart', 1, ...
		'proxfun', f.prox, ...
		'alpha_restart', 0, ...
		'momentum', 0, ...
		'niter', 1000, 'isave', 'all', 'chat', 1);
	xinf = xinf(:,end);
	f.inf = f.cost(xinf);
	%pr log10(f.inf - f.min)
end

%fp = @(x) max(x - f.min, eps^2);
%fp = @(x) x - f.min;
fp = @(x) x - f.inf;
ii = 0:f.niter;
semilogy(ii, fp(cost.gp), '-o', ii, fp(cost.fg), '-+', ...
	ii, fp(cost.ar), '-s', ii, fp(cost.ar90), '-x')
ir_legend({'Gradient Projection', 'FISTA/FPGM', ...
	'FPGM with adaptive restart', 'todo 90'})
