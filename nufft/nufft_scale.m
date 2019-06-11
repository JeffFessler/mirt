 function sn = nufft_scale(Nd, Kd, alpha, beta, Nmid)
%function sn = nufft_scale(Nd, Kd, alpha, beta, Nmid)
%|
%| Compute scaling factors for NUFFT
%|
%| in
%|	Nd,Kd
%|	alpha	{d}
%|	beta	{d}
%|
%| option
%|	Nmid	[d]		midpoint: floor(Nd/2) or default (Nd-1)/2
%|
%| out
%|	sn	[[Nd]]		scaling factors
%|
%| Copyright 2004-7-8, Jeff Fessler, University of Michigan

if nargin == 1 && streq(Nd, 'test'), nufft_scale_test, return, end
if nargin < 4, help(mfilename), help(mfilename), error(mfilename), end

if nargin < 5, Nmid = (Nd-1)/2; end

dd = length(Nd);
if dd == 1 && ~iscell(alpha) % 1D case
	sn = nufft_scale1(Nd(1), Kd(1), alpha, beta, Nmid(1));
return
end


% scaling factors: "outer product" of 1D vectors
sn = 1;
for id=1:dd
	tmp = nufft_scale1(Nd(id), Kd(id), alpha{id}, beta{id}, Nmid(id));
	sn = sn(:) * tmp';
end
if length(Nd) > 1
	sn = reshape(sn, Nd);	% [(Nd)]
else
	sn = sn(:);	% [(Nd)]
end


% Compute scaling factors for 1D NUFFT (from Fourier series coefficients)
% in:
%	N,K
%	alpha
%	beta
% out:
%	sn	[N]		scaling factors
%
% Copyright 2001-10-4, Jeff Fessler, The University of Michigan

function sn = nufft_scale1(N, K, alpha, beta, Nmid)

if ~isreal(alpha(1)), error 'need real alpha_0', end
L = length(alpha) - 1;

%
% compute scaling factors from Fourier coefficients
%
if L > 0
	sn = zeros(N,1);
	n = [0:(N-1)]';
	i_gam_n_n0 = 1i * (2*pi/K) * (n - Nmid) * beta;

	for l1=-L:L
		alf = alpha(abs(l1)+1);
		if l1 < 0, alf = conj(alf); end
		sn = sn + alf * exp(i_gam_n_n0 * l1);
	end

else
	sn = alpha * ones(N,1);
end


% self test
function nufft_scale_test

N = 100;
K = 2*N;
alpha = [1.0 -0.0 -0.2];
sn = nufft_scale(N, K, alpha, 1);
if im
	n = [0:(N-1)]';
	clf, plot(n, real(sn), 'r-', n, imag(sn), 'b-')
	ir_legend({'$s[n]$ real', '$s[n]$ imag'})
	ylabelf('$s[n]$')
	xlabelf('$n$')
end

pr minmax(real(sn))
pr minmax(imag(sn))
