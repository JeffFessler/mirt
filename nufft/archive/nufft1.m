 function [X, times] = nufft1(x, st)
%function [X, times] = nufft1(x, st)
%	SUPERCEDED BY nufft.m
%	Compute NUFFT of 1D signal x
%	in:	x	[N1,L]		set of 1D input signals
%		st			structure precomputed by nufft1_init()
%					st.vm is [K,M]
%	out:	X	[M,L]		set of (approximate) DTFT values
%		times	[2,1]		etime spent in 1: fft, 2: all else.
%
%	Copyright 2001-9-16	Jeff Fessler	The University of Michigan

%
%	if no arguments, then run a simple test
%
if nargin < 1
	help(mfilename)

	N1 = 2^6;
	J1 = 9;
	K1 = 2*N1;
	x = [[1:N1/2]'; ones(N1/2,1)];	% test signal
	rng(0)
	x = randn(N1,1);
	n_shift = 2.7;			% try a random shift just to stress it

	gam = 2*pi/K1;
	o1 = linspace(0,3*gam,N1)';
	printf('err unif %e', ...
		max(nufft1_err_mm(o1, N1, J1, K1, 'sinc', 'uniform')))
	printf('err best %e', ...
		max(nufft1_err_mm(o1, N1, J1, K1, 'sinc', 'best')))

	Xd = dtft1(x, o1, n_shift);
	tic
	st = nufft1_init(o1, N1, J1, K1, n_shift, 'uniform');
	printf('initialization time %g', toc)
	[Xn, times] = nufft1(x, st);
	if 1	% Kaiser-Bessel
		sk = nufft1_init_kb(o1, N1, J1, K1, n_shift, 'kaiser');
		Xk = nufft1(x, sk);
		printf('kb max %% difference = %g', max_percent_diff(Xd,Xk))
	end
	if 1
		sb = nufft1_init(o1, N1, J1, K1, n_shift, 'best');
		Xb = nufft1(x, sb);
		printf('best max %% difference = %g', max_percent_diff(Xd,Xb))
	end
	if 1	% min-max with KB scaling
		smmkb = nufft1_init(o1, N1, J1, K1, n_shift, 'kb');
		Xmmkb = nufft1(x, smmkb);
		printf('mmkb max %% difference = %g', max_percent_diff(Xd,Xmmkb))
	end
	printf('unif max %% difference = %g', max_percent_diff(Xd,Xn))
	clf, plot(	o1/gam, real(Xd), 'y-', o1/gam, imag(Xd), 'g-', ...
			o1/gam, real(Xk), 'm:', o1/gam, imag(Xk), 'r:', ...
			o1/gam, real(Xn), 'b.', o1/gam, imag(Xn), 'c.')
	legend('DTFT r', 'DTFT i', ...
		'KB r', 'KB i', ...
		'NUFFT r', 'NUFFT i')
	printf('time fft=%g other=%g', times(1), times(2))
return
end

time0 = cputime;

if size(x,1) ~= st.N1, error length, end

%
%	apply scaling factors before FFT
%
if size(x,2) > 1
	x = x .* st.sn(:,ones(1,size(x,2)));
else
	x = x .* st.sn;
end

time1 = cputime;

Yk = fft(x, st.K1)';	% [K,L]' oversampled FFT with zero padding at end

time2 = cputime;

X = (Yk * st.vm)';	% [M,L] interpolate using precomputed sparse matrix
%X = st.vm' * Yk;	% this more obvious way is slower!

time3 = cputime;
times(1) = time2-time1;
times(2) = time1-time0 + time3-time2;
