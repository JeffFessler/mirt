 function st = nufft1_init(om, N1, J1, K1, n_shift, alpha, beta)
%function st = nufft1_init(om, N1, J1, K1, n_shift, alpha, beta)
%	Initialize structure for 1D NUFFT.
%	In particular, create the (sparse) interpolation matrix.
%	This interpolator is designed for ordinary FFT (with k=0,...,K-1).
%	in:
%		om	[M,1]	digital frequency omega in radians
%		N1		signal length
%		J1		# of neighbors used per frequency location
%		K1		FFT size (should be > N1)
%		n_shift		n = 0-n_shift to N-1-n_shift.  default is 0.
%		alpha	[L,1]	Fourier series coefficients of scaling factors
%				use 'best' for a good choice
%				or 'kb' for an even better KB-based choice
%		beta		scale gamma=2pi/K by this in Fourier series
%				typically is K/N (me) or 0.5 (Liu)
%	out:
%		st.vm	[K,M]	sparse interpolation matrix
%		st.sn	[N,1]	scaling factors (computed from alpha)
%		st.om,N1,J1,K1	copies of inputs
%
%	Like fft(), the NUFFT expects the signals to be x(0), ..., x(N-1).
%	Use n_shift = N/2 for x(-N/2), ..., x(N/2-1).
%
%	Copyright 2001-8-16	Jeff Fessler	The University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

if ~isvar('n_shift') || isempty(n_shift)
	n_shift = [0];		% default is no shift, so n=0,1,...,N-1
end
if ~isvar('alpha') || isempty(alpha)
	alpha = [1];	% default Fourier series coefficients of scaling factors
end
if ~isvar('beta') || isempty(beta)
	beta = 0.5;	% default is Liu version for now
end

%
%	see if 'best' or 'kb' alpha is desired
%
if ischar(alpha)
	if streq(alpha, 'uniform')
		alpha = 1; beta = 0.5;
	elseif streq(alpha, 'best')
		[alpha, beta, ok] = nufft_best_alpha(J1, 0, K1/N1);
		if ~ok, error 'unknown J,K/N', end
	elseif streq(alpha, 'kb')	% recommended KB-based scaling factors!
		[alpha, beta] = nufft_alpha_kb_fit(N1, J1, K1);
	else
		error 'unknown alpha argument'
	end
end

st.om = om;
st.N1 = N1;
st.J1 = J1;
st.K1 = K1;
st.tol = 0;
st.alpha = alpha;
st.beta = beta;	% scale Fourier series frequency by this

M = length(om);
st.M = M;

gam = 2*pi/K1;
phase_scale = i * gam * (N1-1)/2;	% for phase terms

%
%	compute scaling factors from Fourier coefficients
%
st.sn = nufft_scale(N1, K1, alpha, beta);

%
%	[J,J] precomputed matrix
%
T1 = nufft_T(N1, J1, K1, st.tol, alpha, beta);

%
%	[J,M] r vector
%
[r1, arg] = nufft_r(om, N1, J1, K1, alpha, beta);

%
%	standard interpolator linear phase
%
phase = exp(phase_scale * arg);	% [J,M]
u1 = phase .* (T1 * r1);	% [J,M]

%
%	apply phase due to n_shift
%	(account for Hermitian-tranpose of interpolation coefficients)
%
phase1 = exp(i * n_shift * om);		% [M,1] phase due to n_shift
st.vm = u1 .* phase1(:,ones(1,J1))';	% [J,M]

%
%	build [K,M] sparse matrix
%
koff1 = nufft_offset(om, J1, K1); % [M,1] to leftmost nearest neighbor
kk = mod(outer_sum([1:J1]', koff1'), K1) + 1;		% [J,M] {1,...,K}
ii = [1:M]; ii = ii(ones(J1,1),:);			% [J,M]
st.vm = sparse(kk(:), ii(:), st.vm(:), K1, M);		% [K,M]
