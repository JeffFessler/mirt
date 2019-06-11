% wls_grpr_test.m
% test WLS gradient projection method for various preconditionners
% Copyright Dec. 2000, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi') || 1
%	xtrue = [1 2]';
%	a = -0.5;
%	H = [1 a; a 1];
%	b = [-1 1]';

	G = [1 -1; 0 1];
	G = [1 3; 1 1];
	W = 1;
	H = G' * W * G

	b = H * [-1 1]';
	yi = G' \ b;

%	yi = [-2 1]';
	b = G' * W * yi
end

%
% GP
%
if ~isvar('xgp') || 1
	f.niter = 18;
	xinit = [0.5 0]';
	xinit = [0.5 2]';

	D = 0.3 * diag([1 1]);
	D = diag([0.3 0.1]);
	D = 1.5 / norm(H);
	D = diag([1.3 0.3]) / norm(H);
%	D = inv(H);
	xmin = 0;
	xgp = wls_grpr(xinit, G, W, yi, D, xmin, inf, f.niter);
end


if 1
	xmin = H \ b;
	xcon = [0; b(2) / H(2,2)];

	n1 = 51;
	n2 = 53;
	x1 = linspace(-2, 1, n1);
	x2 = linspace(-1, 2, n2);
	[xx1, xx2] = ndgrid(x1, x2);

	x = [xx1(:) xx2(:)]';
	f = sum(x .* (H * x))/2 - b' * x;
	f = reshape(f, n1, n2);
	f = f';

	if im
		clf
		contour(x1, x2, sqrt(f-min(f(:))))
		hold on
		plot([0 0], minmax(x2), 'w-')
		plot(minmax(x1), [0 0], 'w-')
		plot(xmin(1), xmin(2), 'y*')
		plot(xcon(1), xcon(2), 'yx')
		hold off
		grid
	end
end

if im
	hold on
	plot(xgp(1,:), xgp(2,:), 'g-o')
	hold off
	xlabel x1, ylabel x2
end

% ir_savefig 'fig_wls_grpr'
