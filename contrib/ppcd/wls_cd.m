function [xs,phid,wlsden,timeiter,niter] = wls_cd(A, W, yy, xhat, niter, tol)
%	weighted least squares minimized by coordinate descent algorithm
%	Input
%		A	[nn,np]		system matrix
%		W	[nn,nn]		inverse data covariance
%		yy	[nn,1]		noisy data
%		xhat	[np,1] or [nx,ny]	initial estimate
%		tol			stopping criterion for convergence
%		niter			max # of iterations
%
%	Output
%		xs	[np,niter+1]	estimates each iteration
%		phid	[niter+1,1]	objective decrease
%		wlsden	[np,1]		curvature of surrogate
%		timeiter [niter+1,1]	CPU time each iteration
%		niter			# total iterations
%
% Saowapak Sotthivirat
% 8/21/00

[nx,ny]=size(xhat);
np = ncol(A);

yy = yy(:);
xhat = xhat(:);
xs(:,1) = xhat;

% Precompute the curvature
wlsden = diag(A'*W*A);

wii = diag(W);
res = A * xhat - yy;
phi0 = res'*W*res/2;
ii = 1;
phid(ii) = 0;
t0 = cputime;

while phid(ii) < tol & ii <= niter

	for j=1:np
		xjold = xhat(j);
		xjnew = xjold;
		phijdot = sum(A(:,j).*wii.*res);
		diff = phijdot/wlsden(j);
		xjnew = xjnew - diff;
		xhat(j) = max(0,xjnew);
		res = res - A(:,j) * diff;
	end

	ii = ii+1;
	xs(:,ii) = xhat;
	phid(ii) = phi0 - res'*W*res/2;
	timeiter(ii) = cputime - t0;
end
niter = ii;
