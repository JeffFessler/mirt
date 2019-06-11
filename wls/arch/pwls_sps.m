 function xs = pwls_sps(x, Gt, wi, yi, P, niter, pixmax, ldenom, gi)
%function xs = pwls_sps(x, Gt, wi, yi, P, niter, pixmax, ldenom, gi)
%
%	penalized weighted least squares image reconstruction
%	using separable paraboloidal surrogates algorithm
%	cost function: J(x) = (y-Gx) W (y-Gx) / 2 + x' C' C x / 2
%	Inputs required
%		x	[np,1]		initial estimate
%		Gt	[nd,np]		system matrix', aij >= 0 required!
%		wi	[nd,1]		weighting sinogram
%		yi	[nd,1]		sinogram (of noisy line integrals)
%	Inputs optional
%		P	P.C [nc,np]	penalty matrix
%			P.wpot		penalty weighting function
%		niter			# of iterations
%		ldenom	[np,1]		"likelihood" denominator (optional)
%	Output
%		xs	[np,niter]	iterates
%
%	Copyright Nov. 2000,	Jeff Fessler, University of Michigan

warning 'obsolete.  use pwls_sps_os'

if ~isvar('gi') || isempty(gi)
	gi = sum(Gt)';		% g_i = sum_j g_ij
end

if ~isvar('wi') || isempty(wi)
	wi = ones(size(yi));
end

if ~isvar('pixmax'), pixmax = []; end

%
%	denominator array for likelihood terms
%
if ~isvar('ldenom') || isempty(ldenom)
	ldenom = Gt * (gi .* wi);
end

C = P.C;
%n_per_k = sum(C' ~= 0)';
%if max(n_per_k) == 2
%	warning 'using depierro = 2 version for consistency with ASPIRE'
%	n_per_k(:) = 2;
%end
%	Cfac = (C .* C)' * (n_per_k .* P.wpot(0));	% max curvature
CCnt = (spdiag(sum(C' ~= 0)) * (C .^2))';

%
%	loop over iterations
%
xs = zeros(length(x), niter);
x = bound(x,pixmax);
xs(:,1) = x;
for it=2:niter
	printm('PWLS-SPS iteration %d', it-1)

	li = Gt' * x;
	lgrad = Gt * (wi .* (yi - li));

	Cx = C * x;
	wx = P.wpot(P.wt, Cx);
	pgrad = C' * (wx .* Cx);

	num = lgrad - pgrad;

%	Cfac = (C .* C)' * (n_per_k .* P.wpot(Cx));
	Cfac = CCnt * wx;

	denom = ldenom + Cfac;

	x = x + num ./ denom;	% the update!
	x = bound(x, pixmax);

	xs(:,it) = x;
end

%
%	enforce upper (and possibly lower) bound constraints
%
function z = bound(x, pixmax_min)
if isempty(pixmax_min)
	z = x;
elseif length(pixmax_min) == 1
	z = min(x,pixmax_min);
elseif length(pixmax_min) == 2
	z = min(x,pixmax_min(2));
	z = max(x,pixmax_min(1));
else
	error 'bad pixmax'
end
