 function pgd_step_test
%function pgd_step_test
% test the PGD algorithm using 2D LS cost function

% cost function terms
kap = 4;
A = [1 0; 0 sqrt(kap)];
W = eye(2);
M = eye(2);
yy = 0;

f.niter = 10;
x = [-kap; 1];

% run PSD
xpsd = qpwls_psd(x, A, W, yy, 0, 'precon', 1, 'niter', f.niter, 'isave', 'all');

% run PGD
data = {yy, A, W};
f.step = 3.3; % best step is 0.4, so illustrate suboptimal step
[xpgd1 steps] = pgd_step(x, data, @costgrad, 'precon', M, ...
		'step_method', 'user', 'step', f.step, ...
		'niter', f.niter, 'isave', 'all', 'chat', 0);
	
%pr steps
[xpgd2 steps] = pgd_step(x, data, @costgrad, 'precon', M, ...
		'step_method', 'der2', ...
		'niter', f.niter, 'isave', 'all', 'chat', 0);
%pr steps

if im
	clf, subplot(211)
	plot(	xpsd(1,:), xpsd(2,:), 'y-o', ...
		xpgd1(1,:), xpgd1(2,:), 'b--+', ...
		xpgd2(1,:), xpgd2(2,:), 'm:.')
	xlabel x1, ylabel x2
	title 'QPWLS example'
	leg = {'PSD', sprintf('PGD:%g', f.step), 'PGD:Newton'};
	legend(leg{:}, 3)
end

% cost function map
if ~isvar('qq')
	x1 = linspace(-4.5,0.5,41)';
	x2 = linspace(-1,1.5,43)';
	[xx1 xx2] = ndgrid(x1,x2);
	qq = 0 * xx1;
	for i1=1:length(x1)
		for i2=1:length(x2)
			x = [x1(i1) x2(i2)]';
			qq(i1,i2) = norm(sqrtm(W) * (yy - A * x)).^2;
		end
	end
end

if im
	hold on
	contour(x1, x2, qq', 8)
	plot(0,0, 'rx')
	hold off
	axis equal
	axis([-4.5 0.5 -1 1.5])
	xtick(-4:0)
	ytick(-1:1)
	colormap(0.3+0.7*gray)
end

if im
	subplot(212)
	plot(	0:f.niter, sum(abs(xpsd), 1), 'y-o', ...
		0:f.niter, sum(abs(xpgd1)), 'b--x', ...
		0:f.niter, sum(abs(xpgd2)), 'm:.')
	xlabel 'iteration', ylabel '||x||_1'
	legend(leg{:}, 1)
end


% cost function and gradient for basic WLS
function [cost, grad] = costgrad(x, data)
y = data{1};
A = data{2};
W = data{3};

p = A * x;
cost = 0.5 * norm(sqrtm(W) * (y - p)).^2;
if nargout > 1
	grad = -A' * (W * (y - p));
end
