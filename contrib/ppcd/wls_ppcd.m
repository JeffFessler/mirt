function [xs,phid,wlsden,timeiter,niter] = wls_ppcd(A, W, yy, xhat, nset_col,nset_row, niter, tol)
%	weighted least squares minimized by PPCD algorithm
%
%	Input
%		A	[nn,np]		system matrix
%		W	[nn,nn]		inverse data covariance
%		yy	[nn,1]		noisy data
%		xhat	[nx,ny]		initial estimate
%		nset_col		# sets in column (1 or even)
%		nset_row		# sets in row (1 or even)
%					(total #sets = nset_col*nset_row)
%		tol			stopping criterionfor convergence
%		niter			maximum # of iterations
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

[nx,ny]=size(xhat);
[nn,np] = size(A);

yy = yy(:);
xhat = xhat(:);
xs(:,1) = xhat;

%
%	build sets
%
nset = nset_col*nset_row;
nc = nx/nset_col;
nr = ny/nset_row;
setind	= zeros(np, nset);

xb1 = zeros(ny,nx);
xb1(1:nr,1:nc) = 1;
block1 = find(xb1);
lbk = length(block1);

ind = nx*nc*ones(nset_row,1)*([1:nset_col]-1) + ... \
		[1:nr:nx]'*ones(1,nset_col);
ind = ind(:);
lind = length(ind);

% Each col of Jk corresponds to indexed of pixels in each set.
Jk = block1*ones(1,lind) + ones(lbk,1)*ind' - 1;

% Elements in set1 = 1 if they belong to that set (each col).
set1 = zeros(np,nset);
for k = 1:nset
	set1(Jk(:,k),k) = ones(size(Jk(:,k)));
end

% Precompute rhoik (additive convexity constants) &
% quadratic part of denominator
invrhoik = zeros(np,nset);
wlsden	= zeros(np, 1);
for k = 1:nset
	setk = set1(:,k);
	Jk1 = Jk(:,k);

	num = A*ones(np,1);
	denom = A*setk;
	i = find(denom ~=0);
	j = find(i<=np);
	invrhoik(i(j),k) = num(i(j))./denom(i(j));
	t = diag(A'*W*A).*invrhoik(:,k);
	wlsden(Jk1) = t(Jk1);
end

% compute tik
tik = zeros(nn,nset);
for k=1:nset
	Jk1 = Jk(:,k);
	tik(:,k) = A(:,Jk1)*xhat(Jk1);
end

onek = ones(nset,1);

wik = W*[invrhoik;zeros(nn-np,nset)];

res = A * xhat - yy;
phi0 = res'*W*res/2;
ii = 1;
phid(ii) = 0;
t0 = cputime;

while phid(ii) < tol & ii <= niter

	psidot = W*res;

	for k=1:nset
		wikk = wik(:,k);
		indr = find(wikk);
		qdot = psidot;
		jkk = Jk(:,k);

		for jj=1:lbk

			j = jkk(jj);
			xjold = xhat(j);
			xjnew = xjold;
			Qjdot = A(:,j)' * qdot;
			diff = Qjdot / wlsden(j);
			xjnew = xjnew - diff;
			qdot = qdot - A(:,j).*wikk*diff;
			xhat(j) = max(0,xjnew);
		end
		tik(indr,k) = tik(indr,k) - ...
			(psidot(indr)-qdot(indr))./wikk(indr);
	end

	tk = tik*onek;
	res = tk - yy;
	ii = ii+1;
	xs(:,ii) = xhat;
	phid(ii) = phi0 - res'*W*res/2;
	timeiter(ii) = (cputime - t0)/nset;
end

niter = ii;
