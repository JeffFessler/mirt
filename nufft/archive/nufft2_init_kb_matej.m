function st = nufft2_init_kb(om, N1, N2, J1, J2, K1, K2, n_shift, kernel, useloop, kb_m12, kb_alpha1, kb_alpha2)
%function st = nufft2_init_kb(om, N1, N2, J1, J2, K1, K2, n_shift, kernel, useloop, kb_m12, kb_alpha1, kb_alpha2)
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
%               kb_m12          K-B order parameter (in both directions)
%               kb_alpha1, kb_alpha2  K-B alpha parameter (in each direction)
%	out
%		st.p		[M, N1*N2]	sparse interpolation matrix
%		st.sn		[N1,N2]		scaling factors
%		st.N?,J?,K?,om	copies of inputs
%
%	Like fft(), the NUFFT expects the signals to be x(0,0), ...
%	Use n_shift = [N1/2, N2/2] for x(-N1/2,-N2/2), ...
%
%	Copyright 2002-4-12	Jeff Fessler	The University of Michigan

kb_m1 = kb_m12;	% for simplicity
kb_m2 = kb_m12;

kernel1 = kaiser_bessel('inline', J1, kb_alpha1, kb_m1);
kernel2 = kaiser_bessel('inline', J2, kb_alpha2, kb_m2);

% the following form will use the numerical FT as scaling factors:
st = nufft_init(om, [N1 N2], [J1 J2], [K1 K2], n_shift, ...
	{kernel1, kernel2})

% the following form will use the analytical FT as scaling factors:
st = nufft_init(om, [N1 N2], [J1 J2], [K1 K2], n_shift, ...
	'kaiser', [kb_alpha1 kb_alpha2], [kb_m1 kb_m2])
