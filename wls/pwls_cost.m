 function [cost, fit, reg] = pwls_cost(xs, A, W, yi, R, mask)
%function [cost, fit, reg] = pwls_cost(xs, A, W, yi, R, mask)
%|
%| compute PWLS cost for each column of xs
%|
%| in
%|	xs	[np niter]	iterates
%|	A	[nd np]		system matrix
%|	W	[nd nd]		data weighting matrix, usually diag_sp(wi)
%|	yi	[nd 1]		data
%|	R			penalty object (see Reg1.m or Robject.m)
%|				with penalty cost method: R.penal(R, x)
%|				or just *sparse* C matrix (quadratic only)
%|	mask	[nx ny ...]	optional mask, iff xs is [nx ny ... niter]
%|
%| out
%|	cost	[niter 1]	cost
%|	fit	[niter 1]	(y-A*x)'W(y-Ax)/2
%|	reg	[niter 1]	R(x)
%|
%| Copyright 2002-2-12, Jeff Fessler, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

if ~isvar('R'), R = []; end

if isvar('mask') && ~isempty(mask)
	xs = reshapee(xs, numel(mask), []);	% [(*N) niter]
	xs = xs(mask(:), :);			% [np niter]
end

niter = size(xs,2);
reg = zeros(niter,1);

if isempty(R)
	warning 'empty R means no penalty'

elseif issparse(R) || isa(R, 'Fatrix')
	C = R; % trick!
	if size(C,2) == size(C,1)
		warning 'square C is quite unusual!?'
	end
	for ii=1:niter
		reg(ii) = sum(abs(C * xs(:,ii)).^2)/2;
	end

elseif isstruct(R) || isa(R, 'strum')
	for ii=1:niter
		reg(ii) = R.penal(R, xs(:,ii));
	end

else
	keyboard
	error 'bad R'
end

fit = zeros(niter,1);
for ii=1:niter
	resid = yi - A * xs(:,ii); % predicted measurements
	fit(ii) = resid' * (W * resid) / 2;
end

cost = fit + reg;
cost = reale(cost, 'warn'); % trick: x'*W*x is not always real for complex values

if ~nargout
	pr fit
	pr reg
	pr cost
end
