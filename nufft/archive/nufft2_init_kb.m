 function st = nufft2_init_kb(om, N1, N2, J1, J2, K1, K2, n_shift, kernel, useloop)
%function st = nufft2_init_kb(om, N1, N2, J1, J2, K1, K2, n_shift, kernel, useloop)
%	Initialize structure for 2D NUFFT using KB interpolator,
%	particularly the interpolation matrix in sparse format.
%	in
%		om [M,2]	frequencies in radians
%		N1,N2		image dimensions
%		J1,J2		# of neighbors used (in each direction)
%		K1,K2		FFT sizes (should be > N1,N2)
%		n_shift [2]	n = 0-n_shift to N-1-n_shift
%		kernel		inline kernel function, or 'kaiser'
%		useloop		1 for slow version with less memory 
%	out
%		st.p		[M, N1*N2]	sparse interpolation matrix
%		st.sn		[N1,N2]		scaling factors
%		st.N?,J?,K?,om	copies of inputs
%
%	Like fft(), the NUFFT expects the signals to be x(0,0), ...
%	Use n_shift = [N1/2, N2/2] for x(-N1/2,-N2/2), ...
%
%	Copyright 2002-4-12	Jeff Fessler	The University of Michigan

if nargin < 7, help(mfilename), error args, end
if ~isvar('n_shift') || isempty(n_shift), n_shift = [0 0]; end
if ~isvar('kernel') || isempty(kernel), kernel = 'kaiser'; end
if ~isvar('useloop') || isempty(useloop), useloop = false; end

%
%	see if 'best' alpha is desired
%
is_kaiser = false;
if ischar(kernel)
	if streq(kernel, 'kaiser')
		is_kaiser = true;
		[st.kernel1 st.kb_alf1 st.kb_m1] = kaiser_bessel('inline', J1);
		[st.kernel2 st.kb_alf2 st.kb_m2] = kaiser_bessel('inline', J2);
	else
		error 'unknown kernel'
	end
elseif isa(kernel, 'inline')
	st.kernel1 = kernel;
	st.kernel2 = kernel;
else
	error 'kernel must be string or inline'
end

st.tol	= 0;

st.J1 = J1;
st.J2 = J2;
st.N1 = N1;
st.K1 = K1;
st.N2 = N2;
st.K2 = K2;

st.om = om;
M = size(om,1);
st.M = M;
if size(om,2) ~= 2, error 'omega needs two columns', end

%
%	scaling factors
%
if is_kaiser
	n1 = [0:N1-1]'-(N1-1)/2;
	n2 = [0:N2-1]'-(N2-1)/2;
	st.sn1 = 1 ./ kaiser_bessel_ft(n1/K1, J1, st.kb_alf1, st.kb_m1, 1);
	st.sn2 = 1 ./ kaiser_bessel_ft(n2/K2, J2, st.kb_alf2, st.kb_m2, 1);
else
	st.sn1 = 1 ./ nufft_interp_zn(0, N1, J1, K1, st.kernel1);
	st.sn2 = 1 ./ nufft_interp_zn(0, N2, J2, K2, st.kernel2);
%else
%	st.sn1 = nufft_scale(N1, K1, alpha1, beta1);
%	st.sn2 = nufft_scale(N2, K2, alpha2, beta2);
end
st.sn = st.sn1 * st.sn2';	% outer product

%
%	[J,M] coefficient vector
%
[c1, arg1] = nufft_coef(om(:,1), J1, K1, st.kernel1);	% [J1,M]
[c2, arg2] = nufft_coef(om(:,2), J2, K2, st.kernel2);	% [J2,M]

%
%	1D interpolation coefficient vectors.  will need kron of these later.
%
gam1 = 2*pi/K1;
gam2 = 2*pi/K2;
phase_scale1 = i * gam1 * (N1-1)/2;
phase_scale2 = i * gam2 * (N2-1)/2;

phase1 = exp(phase_scale1 * arg1);	% [J1,M] linear phase
phase2 = exp(phase_scale2 * arg2);	% [J2,M] linear phase
u1 = phase1 .* c1;			% [J1,M]
u2 = phase2 .* c2;			% [J2,M]
clear c1 c2 arg1 arg2 phase1 phase2 phase_scale1 phase_scale2

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
