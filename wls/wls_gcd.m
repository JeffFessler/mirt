  function [xs, wlsden] = wls_gcd(A, W, yy, xhat, niter, ng, chat)
%|function [xs, wlsden] = wls_gcd(A, W, yy, xhat, niter, ng, chat)
%| weighted least squares minimized by
%| grouped coordinate descent algorithm
%| in
%|	A	[nn np]			system matrix
%|	W	[nn nn]			inverse data covariance
%|	yy	[nn 1]			noisy data
%|	nder1	[nn 1]
%|	xhat	[np 1] or [nx ny]	initial estimate
%|	niter		# total iterations
%|	ng		# groups, or, if array, groups
%| out
%|	xs	[np niter]		estimates each iteration
%|
%| Copyright July 1996, Jeff Fessler, University of Michigan

if nargin < 6, help(mfilename), error(' '), end

if 1 ~= exist('chat'), chat = 1; end
%chat = (1==exist('mask'));

[nx,ny] = size(xhat);

if 1~=exist('mask') || isequal(mask, 1)
	mask = ones(nx,ny);
end

np = ncol(A);

if numel(xhat) ~= np
	xhat = xhat(mask(:));
end

if ~isempty(nder1)
	back_nder1 = A' * nder1(:);
else
	back_nder1 = zeros(np,1);
end
yy = yy(:);
xs = zeros(np, niter);
xs(:,1) = xhat;

if chat
	im(321, embed(xhat,mask), 'x0')
end

% build groups
if numel(ng) == 1
	groups	= zeros(np, ng);
	for ig = 1:ng
		gg = [ig:ng:np]';
		groups(gg,ig) = ones(size(gg));
	end
elseif numel(ng) == 2
	groups	= group2d(nx, ny, ng(1), ng(2));
	ng	= ncol(groups);
else
	groups	= ng;
	ng	= ncol(groups);
end

% overlapping groups?
if any(sum(groups')) > 1, error 'overlap groups', end

if chat
	im(322, embed(groups * [1:ng]',mask), 'Groups')
end

% precompute quadratic part of denominator
% change '1' to 'ig' for overlapping groups
wlsden	= zeros(np, 1);
if 1
	for ig = 1:ng
		gg = groups(:,ig);
		printm('precompute %g, assuming nonnegative A', ig)

%		wlsden(gg,1) = abs(A(:,gg))' * (W * sum(abs(A(:,gg)'),1)');
%		wlsden(gg,1) = A(:,gg)' * (W * sum(A(:,gg)')',1);

		t = A * gg;	% faster than: sum(A(:,gg)',1)';
		t = W * t;	% only slightly slower than: diag(W) .* t;
		%	faster than either: A(:,gg)'*t or A' * t !
		t = t' * A;
		wlsden(gg,1) = t(:,gg)';
	end
else
	wlsden	= ones(np, 1);	% for testing
end

if chat
	im(323, embed(wlsden,mask), 'WLS Denominator')
end

% compare GCA and SCA denominators (for WLS part)
if chat && (np < 1000)
	t = full(sum(A .* (W * A))'); % SCA den
	t = t(wlsden ~= 0) ./ wlsden(wlsden ~= 0);
	if min(t(:)) < 1-10*eps
		printm(['Under-relax. range:' sprintf(' %g', minmax(t))])
		t = embed(t,wlsden);
		im(324, embed(t,mask), 'Under-relaxer')
	else
		printm 'Note: GCA with ideal denominator!'
	end
end

% 'depierro' factors for penalty? (GCA vs SCA)
for ig = 1:ng
	gg = groups(:,ig);
	t = abs(C(:,gg))' * sum(abs(C(:,gg)'),1)';
	t = sum(C(:,gg).^2)' ./ t; % SCA / GCA
	t = minmax(t);
	if min(t) < 1-10*eps
		printm(['Virtual under-relax. range:' sprintf(' %g', t)])
		error 'depierro not done'
	end
	clear t
end

if chat
	im(325, embed(sum(C.^2)',mask), 'sum(C.^2)')
	input('hit key');
end

res = A * xhat - yy;

for ii = 2:niter

	for ig = 1:ng
		gg = groups(:,ig);

		xn = xhat(gg);
		xw = xn;

		ader = ( (W * res)' * A(:,gg) )' + back_nder1(gg);

		for jj = 1:nsub
			dx = C * xhat;
			eval(sprintf(wstring, 'wt', 'dx'));
			denom = wlsden(gg) + (C(:,gg).^2)' * wt;
			num = ader + wlsden(gg).*(xw-xn) ...
				+ C(:,gg)' * (wt .* dx);
			xw = xw - num ./ denom;
			xhat(gg) = xw;
		end

		res = res + A(:,gg) * (xw - xn);
	end

	xs(:,ii) = xhat;
	if chat && (ii < 10 || rem(ii,10) == 0)
		im(326, embed(xhat,mask), 'x hat')
		if ny == 1
			im(324, wt, 'wts')
		end
	end
	if chat
		printm('max %g %g', max(xhat), max(abs(xhat-xs(:,ii-1)))/max(xhat)*100)
	end
end
