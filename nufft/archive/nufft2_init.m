 function st = nufft2_init(om, N1, N2, J1, J2, K1, K2, n_shift, useloop, alpha, beta)
%function st = nufft2_init(om, N1, N2, J1, J2, K1, K2, n_shift, useloop, alpha, beta)
%	initialize structure for 2D NUFFT
%	particularly the interpolation matrix in sparse format
%	in:
%		om [M,2]	frequencies in radians
%		N1,N2		image dimensions
%		J1,J2		# of neighbors used (in each direction)
%		K1,K2		FFT sizes (should be > N1,N2)
%		n_shift [2]	n = 0-n_shift to N-1-n_shift
%		useloop		set to 1 for slow, low memory, loop
%		alpha		scaling factors
%				use 'kb' for very good KB-based scaling factors
%		beta		frequency parameter for scaling factors
%	out:
%		st.p		[M, N1*N2]	sparse interpolation matrix
%		st.sn		[N1,N2]		scaling factors
%		st.N1,J1,K1	copies of inputs
%		st.N2,J2,K2	""
%		st.om		""
%
%	Like fft(), the NUFFT expects the signals to be x(0,0), ...
%	Use n_shift = [N1/2, N2/2] for x(-N1/2,-N2/2), ...
%
%	Copyright 2001-9-17	Jeff Fessler	The University of Michigan

if nargin < 7, help(mfilename), error args, end
if ~isvar('n_shift') || isempty(n_shift), n_shift = [0 0]; end
if ~isvar('useloop') || isempty(useloop), useloop = false; end
if ~isvar('alpha') || isempty(alpha), alpha = 1; end
if ~isvar('beta') || isempty(beta), beta = 0.5; end

alpha1 = alpha;
alpha2 = alpha;
beta1 = beta;
beta2 = beta;

%
%	see if 'best' or 'kb' alpha is desired
%
if ischar(alpha)
	if streq(alpha, 'best')
		[alpha1, beta1, ok] = nufft_best_alpha(J1, 0, K1/N1);
		if ~ok, error 'unknown J,K/N', end
	elseif streq(alpha, 'kb')
		[alpha1, beta1] = nufft_alpha_kb_fit(N1, J1, K1);
	else
		error 'unknown alpha argument'
	end

	if streq(alpha, 'best')
		[alpha2, beta2, ok] = nufft_best_alpha(J2, 0, K2/N2);
		if ~ok, error 'unknown J,K/N', end
	elseif streq(alpha, 'kb')
		[alpha2, beta2] = nufft_alpha_kb_fit(N2, J2, K2);
	else
		error 'unknown alpha argument'
	end
end
clear alpha beta

st.alpha1 = alpha1;
st.alpha2 = alpha2;
st.beta1 = beta1;
st.beta2 = beta2;
st.tol	= 0;

st.N1 = N1;
st.J1 = J1;
st.K1 = K1;

st.N2 = N2;
st.J2 = J2;
st.K2 = K2;

st.om = om;
M = size(om,1);
st.M = M;
if size(om,2) ~= 2, error 'omega needs two columns', end

%
%	scaling factors
%
st.sn1 = nufft_scale(N1, K1, alpha1, beta1);
st.sn2 = nufft_scale(N2, K2, alpha2, beta2);
st.sn = st.sn1 * st.sn2';	% outer product

%
%	[J,J] precomputed matrix
%
T1 = nufft_T(N1, J1, K1, st.tol, alpha1, beta1);	% [J1,J1]
T2 = nufft_T(N2, J2, K2, st.tol, alpha2, beta2);	% [J2,J2]

[r1, arg1] = nufft_r(om(:,1), N1, J1, K1, alpha1, beta1);	% [J1,M]
[r2, arg2] = nufft_r(om(:,2), N2, J2, K2, alpha2, beta2);	% [J2,M]

%
%	1D interpolation coefficient vectors.  will need kron of these later.
%
gam1 = 2*pi/K1;
gam2 = 2*pi/K2;
u1 = exp((i * gam1 * (N1-1)/2) * arg1);	% [J1,M] linear phase
u2 = exp((i * gam2 * (N2-1)/2) * arg2);	% [J2,M] linear phase
u1 = u1 .* (T1 * r1);			% [J1,M]
u2 = u2 .* (T2 * r2);			% [J2,M]
clear r1 r2 arg1 arg2

%
%	indices into oversampled FFT components
%
koff1 = nufft_offset(om(:,1), J1, K1);	% [M,1] to leftmost near nbr
koff2 = nufft_offset(om(:,2), J2, K2);

k1 = mod(outer_sum([1:J1]', koff1'), K1) + 1;	% [J1,M] {1,...,K1}
k2 = mod(outer_sum([1:J2]', koff2'), K2) + 1;	% [J2,M] {1,...,K2}
clear koff1 koff2

%
%	build sparse matrix that is [M,K1*K2]
%	with J1*J2 nonzero entries per frequency point
%

%
%	slowest but simplest, for development and testing
%
if useloop
	%
	%	P will be [M,K1*K2], but allocate transpose for now
	%	J1*J2 nonzero entries per frequency point
	%
	P = spalloc(K1*K2, M, J1*J2*M);

	%
	%	loop over desired frequency points
	%
	for ii = 1:M
		kk = outer_sum(k1(:,ii), (k2(:,ii)-1) * K1);
		P(kk(:),ii) = kron(u2(:,ii), u1(:,ii));
	end

	%
	%	apply phase shift
	%	pre-do Hermitian transpose of interpolation coefficients
	%
	phase = exp(i * (om * n_shift(:)));
	st.p = spdiag(phase) * P';

%
%	this version avoids the loop but uses O(J1*J2*M) memory
%
else
	ll1 = reshape(k1, [J1 1 M]);		% [J1,1,M] from [J1,M]
	ll1 = ll1(:,ones(J2,1),:);		% [J1,J2,M], emulating ndgrid
	ll2 = reshape(k2, [1 J2 M]);		% [1,J2,M] from [J2,M]
	ll2 = ll2(ones(J1,1),:,:);		% [J1,J2,M], emulating ndgrid
	kk = ll1 + (ll2 - 1) * K1;		% [J1,J2,M]
	clear l1 l2 ll1 ll2

	uu1 = reshape(u1, [J1 1 M]);		% [J1,1,M] from [J1,M]
	uu1 = uu1(:,ones(J2,1),:);		% [J1,J2,M], emulating ndgrid
	uu2 = reshape(u2, [1 J2 M]);		% [1,J2,M] from [J2,M]
	uu2 = uu2(ones(J1,1),:,:);		% [J1,J2,M], emulating ndgrid
	uu = uu1 .* uu2;			% [J1,J2,M]
	clear u1 u2 uu1 uu2

	%
	%	apply phase shift
	%	pre-do Hermitian transpose of interpolation coefficients
	%
	phase = exp(i * (om * n_shift(:))).';	% [1,M]
	uu = reshape(uu, [J1*J2, M]);			% [J1*J2,M]
	uu = conj(uu) .* phase(ones(1,J1*J2),:);	% [J1*J2*M]

	ii = [1:M]; ii = ii(ones(J1*J2,1),:);	% [J1,J2,M]
	st.p = sparse(ii(:), kk(:), uu(:), M, K1*K2);	% sparse matrix
end
