  function xs = qpwls_art2(x, Gt, wi, y, Ct, epsi, niter, gnorms)
%|function xs = qpwls_art2(x, Gt, wi, y, Ct, epsi, niter, gnorms)
%|
%| identity penalized weighted least squares reconstruction
%| using ART2 algorithm from Censor and Zenios, modified somewhat per Dax
%| cost(x) = |y-Gx|_W^2 / 2 + |Cx|^2 / 2 + eps |x|^2 / 2
%|
%| in
%|	x	[np 1]		initial guess
%|	CAUTION: x must be in the range of G' diag(sqrt(wi)) !!!!!
%|	Gt	[np nd]		transpose of system matrix
%|	wi	[nb na]		weights
%|	y	[nb na]		measurement
%|	Ct	[np nc]		transpose of regularization matrix
%|	epsi	[1]		small parameter
%|	niter			# of iterations
%|	gnorms	[nd 1]		see below
%|
%| out
%|	xs	[np niter]	iterates
%|
%| Copyright 2000-06, Jeff Fessler, University of Michigan

if nargin < 5, help(mfilename), error(mfilename), end

if ~isvar('niter'), niter = 1; end

[nb na] = size(y);
starts = subset_start(na);

% For WLS, premultiply y and postmultiply Gt by W^{1/2}
Wh = spdiag(sqrt(wi(:)), 'nowarn');
y = Wh * y(:);
try
Gt = Gt * Wh;
catch % todo
keyboard
end

% weighted row norms
if ~isvar('gnorms')
	gnorms = sum(Gt.^2); % | e_i' G |^2
	cnorms = sum(Ct.^2);
end

%vg = zeros(ncol(Gt),1);	% G component of v vector
vg = -y;			% trick: G component of v vector (-y)
vc = zeros(ncol(Ct),1);		% C component of v vector

xs(:,1) = x;

iclist = 1:length(vc);
iclist = iclist(cnorms(iclist) ~= 0);

iglist = col(outer_sum(1:nb, (starts-1)*nb));
iglist = iglist(gnorms(iglist) ~= 0);

gdenom = gnorms + epsi;
cdenom = cnorms + epsi;

for it=2:niter
	ticker(mfilename, it, niter)

	% system matrix part
	for ii=iglist'
		[j ignore g] = find(Gt(:,ii));
		xj = x(j);
		step = g' * max(xj,0) + vg(ii);
		step = step / gdenom(ii);
		x(j) = xj - step * g;
		vg(ii) = vg(ii) - step * epsi;
	end

	% regularization part
	% fix: what ordering would be logical here?
	for ii=iclist
		[j ignore c] = find(Ct(:,ii));
		xj = x(j);
		step = c' * max(xj,0) + vc(ii);
		step = step / cdenom(ii);
		x(j) = xj - step * c;
		vc(ii) = vc(ii) - step * epsi;
	end

	xs(:,it) = x;
end
