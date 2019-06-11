 function xs = ir_ip_unwrap_sqs(x, yi, R, niter)
%function xs = ir_ip_unwrap_sqs(x, yi, R, niter)
%|
%| separable quadratic surrogates phase unwrapping (x denotes phase)
%| can be 1d, 2d, etc., depending on R.
%|
%| cost(x) = |y|^2 (1 - cos(x - \angle(y)) + R(x)
%|
%| in
%|	x	[np 1]		initial estimate
%|	yi	[np 1]		measurements (noisy sinogram)
%|	R			penalty object (see Reg1.m)
%|	niter			# of iterations
%| out
%|	xs	[np niter]	iterates
%|
%| Copyright 2002-6-15, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

if ~isvar('niter') || isempty(niter), niter = 1; end
if ~isvar('chat') || isempty(chat), chat = false; end

if ~isvar('R') || isempty(R)
	pgrad = 0; % unregularized default
	Rdenom = 0;
end

ang = angle(yi);
mag2 = abs(yi).^2;

% loop over iterations

xs = zeros(numel(x), niter);
xs(:,1) = x;

for iter = 2:niter
	if chat, printf('unwrap iteration %d', iter-1), end

	s = x - ang;
	grad = mag2 .* sin(s);

	if ~isempty(R)
		pgrad = R.cgrad(R, x);
		Rdenom = R.denom(R, x);
	end

	sr = mod(s+pi,2*pi) - pi; % [-pi,pi]
	denom = mag2 .* nufft_sinc(sr / pi); % curvatures

	num = grad + pgrad;
	den = denom + Rdenom;

	x = x - num ./ den; % relaxed update

	if chat, printf('Range %g %g', min(x), max(x)), end
	xs(:,iter) = x;
end
