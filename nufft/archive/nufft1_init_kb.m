 function st = nufft1_init_kb(om, N1, J1, K1, n_shift, kernel, karg)
%function st = nufft1_init_kb(om, N1, J1, K1, n_shift, kernel, karg)
%	Initialize structure for 1D NUFFT using Kaiser-Bessel interpolator.
%	In particular, create the (sparse) interpolation matrix.
%	This interpolator is designed for ordinary FFT (with k=0,...,K-1).
%	in:
%		om	[M,1]	digital frequency omega in radians
%		N1		signal length
%		J1		# of neighbors used per frequency location
%		K1		FFT size (should be > N1)
%		n_shift		n = 0-n_shift to N-1-n_shift.  default is 0.
%		kernel		inline kernel function or 'kaiser'
%		karg		optional kernel arguments
%	out:
%		st.vm	[K,M]	sparse interpolation matrix
%		st.sn	[N,1]	scaling factors (computed from alpha)
%		st.om,N1,J1,K1	copies of inputs
%
%	Like fft(), the NUFFT expects the signals to be x(0), ..., x(N-1).
%	Use n_shift = N/2 for x(-N/2), ..., x(N/2-1).
%
%	Copyright 2002-4-11	Jeff Fessler	The University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

if ~isvar('n_shift') || isempty(n_shift)
	n_shift = [0];		% default is no shift, so n=0,1,...,N-1
end
if ~isvar('kernel') || isempty(kernel)
	kernel = 'kaiser';
end


%
%	various ways of specifying the kernel
%
is_kaiser = false;
if ischar(kernel)
	if streq(kernel, 'kaiser')
		[st.kernel st.alpha st.kb_m] = kaiser_bessel('inline', J1);
		is_kaiser = true;
	else
		error(sprintf('kernel "%s" not done', kernel))
	end

elseif isa(kernel, 'inline')
	tmp = struct(kernel);
	if ~isempty(findstr(tmp.expr), 'kaiser_bessel')
		warning 'using non-FT scaling factors for KB!?'
%		is_kaiser = false;
	end
	st.kernel = kernel;

else
	error 'kernel must be string or inline'
%	st.kernel = kaiser_bessel('inline', J1, alpha, kb_m);
end

st.om = om;
st.N1 = N1;
st.J1 = J1;
st.K1 = K1;
st.tol = 0;

M = length(om);
st.M = M;

gam = 2*pi/K1;
phase_scale = i * gam * (N1-1)/2;	% for phase terms

%
%	compute KB scaling factors from analytical Fourier transform
%
if is_kaiser
	n1 = [0:N1-1]'-(N1-1)/2;
	st.sn = 1 ./ kaiser_bessel_ft(n1/K1, J1, st.alpha, st.kb_m, 1);
else
	st.sn = 1 ./ nufft_interp_zn(0, N1, J1, K1, st.kernel);
end

%
%	[J,M] coef vector
%
[coef, arg] = nufft_coef(om, J1, K1, st.kernel);

%
%	standard interpolator linear phase
%
phase = exp(phase_scale * arg);	% [J,M]
u1 = phase .* coef;		% [J,M]

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
