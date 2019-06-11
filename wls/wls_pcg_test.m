% wls_pcg_test.m
%
% Test WLS-PCG with a toy complex problem
%
% Copyright 2001-07-06, Jeff Fessler, University of Michigan

% generate random toy data: y=A*x+noise
if 1 || ~isvar('yi')
	im clf, rng(0)
	np = 6;	nd = 8;
	xtrue = randn(np,2) * [1; 1i];			% complex truth
	A = randn(nd,np) + 1i * randn(nd,np);		% complex system
	yi = A * xtrue + 0.1 * randn(nd,2) * [1; 1i];	% complex noise
	W = diag_sp([1:nd]'); % weighted
end

% regularization matrix
if 1 || ~isvar('C')
	nc = 2 * np;
	C = randn(nc,np) + 1i * randn(nc,np);
end

% exact solution
if 1 || ~isvar('xhat')
	xhat = [A'*(W*A) + C'*C] \ (A'*(W*yi));
end

% uniform initial image
xinit = zeros(size(xtrue));
f.niter = 10;

% run PCG
if 1 || ~isvar('xpcg')
	[xpcg info] = qpwls_pcg(xinit, A, W, yi, 0, C, 1, f.niter);
end

equivs(xpcg(:,end), xhat) % check for match

% plot cost function
if 1 && im
	cost.pcg = pwls_cost(xpcg, A, W, yi, sparse(C));
	ii = 0:f.niter-1;
	im clf, subplot(211)
	plot(ii, cost.pcg, '-*')
	xlabel 'iteration', ylabel 'cost function'
	legend 'PCG cost'
end

% plot distance to "exact" xhat
if 1 && im
%	nrms = sqrt(2*pwls_cost(xpcg, 1, 1, xhat, 0)) / norm(xhat);
	subplot(212)
	semilogy(ii, nrms(xpcg, xhat), '-.')
	xlabel 'iteration', ylabel 'NRMS error'
end
