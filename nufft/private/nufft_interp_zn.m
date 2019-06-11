 function zn = nufft_interp_zn(alist, N, J, K, func, Nmid)
%function zn = nufft_interp_zn(alist, N, J, K, func, Nmid)
%|
%| Compute the "zn" terms for a conventional "shift-invariant" interpolator
%| as described in T-SP paper.  Needed for error analysis and for user-
%| defined kernels since I don't provide a means to put in an analytical
%| Fourier transform for such kernels.
%|
%| in
%|	alist	[M]	omega / gamma (fractions) in [0,1)
%|	func		func(k,J) support limited to [-J/2,J/2)
%|				interpolator should not include the linear
%|				phase term.  This routine provides it.
%| option
%|	Nmid	[1]	midpoint: floor(N/2) or default: (N-1)/2
%|
%| out
%|	zn	[N M]	reciprocal of scaling factors
%|
%| Copyright 2001-12-11, Jeff Fessler, University of Michigan

%| zn = \sum_{j=-J/2}^{J/2-1} exp(i gam (alf - j) * n) F1(alf - j)
%| = \sum_{j=-J/2}^{J/2-1} exp(i gam (alf - j) * (n-Nmid)) F0(alf - j)

if nargin == 1 && streq(alist, 'test'), nufft_interp_zn_test, return, end
if nargin < 5, ir_usage, end
if nargin < 6
	Nmid = (N-1)/2; % default: old version
end

gam = 2*pi/K;

if any(alist < 0 | alist > 1), warn 'alist exceeds [0,1]', end

% natural phase function. trick: force it to be 2pi periodic
%Pfunc = @(om, N) exp(-1i * mod0(om,2*pi) * (N-1)/2);

if ~rem(J,2) % even
	jlist = [(-J/2+1):J/2]';
else	% odd
	jlist = [-(J-1)/2:(J-1)/2]';
	alist(alist > 0.5) = 1 - alist(alist > 0.5); % force symmetry!
end

n0 = [0:(N-1)]' - Nmid; % include effect of phase shift!
[nn0, jj] = ndgrid(n0, jlist); % [N J]
zn = zeros(N, length(alist));

for ia=1:length(alist)
	alf = alist(ia);
	jarg = alf - jj;			% [N J]
	e = exp(1i * gam * jarg .* nn0);	% [N J]

	F = func(jarg, J);			% [N J]
	zn(:,ia) = sum(F .* e, 2);		% [N]
end


function nufft_interp_zn_test
alist = [0:19]/20;
N = 2^7;
K = 2 * N;

% linear
J = 4;
func = @(k, J) (1 - abs(k/(J/2))) .* (abs(k) < J/2);

z = nufft_interp_zn(alist, N, J, K, func);

% plot interpolator
if im
	k = linspace(-J/2-1, J/2+1, 101);
	clf, jf pl 1 3
	jf sub 1, plot(k, func(k, J))
	xlabel k, ylabel F_0(k), axis tight

	jf sub 2, plot(1:N, real(z)), axis tight, ylabel 'real(z_n)'
	jf sub 3, plot(1:N, imag(z)), axis tight, ylabel 'imag(z_n)'
end
