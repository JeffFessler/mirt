function [xs,phid,wlsden,timeiter,niter] = wls_gcd(A, W, yy, xhat, ng, niter, tol)
%	weighted least squares minimized by
%	grouped coordinate descent algorithm
%	Input
%		A	[nn,np]		system matrix
%		W	[nn,nn]		inverse data covariance
%		yy	[nn,1]		noisy data
%		xhat	[np,1] or [nx,ny]	initial estimate
%		ng			# groups, or, if array, groups
%		tol			stopping criterion for convergence
%		niter			max
%
%	Output
%		xs	[np,niter+1]	estimates each iteration
%		phid	[niter+1,1]	objective decrease
%		wlsden	[np,1]		curvature of surrogate
%		timeiter [niter+1,1]	CPU time each iteration
%		niter			# total iterations
%		
% Saowapak Sotthivirat
% 8/15/00

[nx,ny] = size(xhat);

np = ncol(A);
yy = yy(:);
xhat = xhat(:);
xs(:,1) = xhat;

%
%	build groups
%
if numel(ng) == 1
	groups = zeros(np, ng);
	for ig = 1:ng
		gg(:,ig) = [ig:ng:np]';
		groups(gg(:,ig),ig) = ones(size(gg(:,ig)));
	end
elseif numel(ng) == 2
	groups	= group2d(nx, ny,ng,ones(ng(1),ng(2)),0);
	ng	= ncol(groups);
	for ig = 1:ng
		gg(:,ig) = find(groups(:,ig).*[1:np]');
	end
else
	groups	= ng;
	ng	= ncol(groups);
end


%	precompute quadratic part of denominator
%	change '1' to 'ig' for overlapping groups
wlsden = zeros(np, 1);
for ig = 1:ng
	gr = groups(:,ig);
	ggi = gg(:,ig);

	t = A * gr;	% faster than: sum(A(:,gg)',1)';
	t = W * t;	% only slightly slower than: diag(W) .* t;
			% faster than either: A(:,gg)'*t or A' * t !
	t = t' * A;
	wlsden(ggi) = t(ggi)';
end

wii = diag(W);
res = A * xhat - yy;
phi0 = res'*W*res/2;
ii = 1;
phid(ii) = 0;
timeiter(ii) = 0;
t0 = cputime;

while phid(ii) < tol & ii <= niter
	
	for ig = 1:ng
		ggi = gg(:,ig);
		xold = xhat(ggi);
		xnew = xold;
		ader = ((W*res)'*A(:,ggi))';
		diff = ader./wlsden(ggi);
		xnew = xnew - diff;
		xhat(ggi) = max(0,xnew);
		res = res - A(:,ggi) * diff;

	end
	ii = ii+1;
	xs(:,ii) = xhat;
	phid(ii) = phi0 - res'*W*res/2;
	timeiter(ii) = cputime - t0;
end

niter = ii;
