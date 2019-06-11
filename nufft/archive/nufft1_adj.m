 function x = nufft1_adj(X, st)
%function x = nufft1_adj(X, st)
%	SUPERCEDED BY nufft_adj.m
%	Apply adjoint of 1D NUFFT operator to X
%	in:	X	[M,L]	set of input spectral vectors
%		st		structure precomputed by nufft1_init()
%	out:	x	[N1,L]	result of applying adjoint
%
%	caution: the adjoint is not the inverse, so x is not really x
%
%	Copyright 2001-9-16	Jeff Fessler	The University of Michigan

%
%	if no arguments, then run a simple test
%
if nargin < 1
	help(mfilename)

	N1 = 8;
	o1 = 2*pi * [0.1 0.3 0.4 0.7 0.9]';
	M1 = length(o1);
	n_shift = 3.3;		% random shift to stress it
	st = nufft1_init(o1, N1, 6, 2*N1, n_shift, 'best');
	A = zeros(M1, N1);
	for jj=1:N1
		x = zeros(N1,1); x(jj,1) = 1;
		A(:,jj) = nufft1(x, st);
	end

	Aa = zeros(N1, M1);
	for ii=1:M1
		y = zeros(M1,1); y(ii,1) = 1;
		Aa(:,ii) = nufft1_adj(y, st);
	end

	printf('max %% difference = %g', max_percent_diff(A,Aa'))
	clear x
return
end


if any(size(X,1) ~= st.M), error length, end

% adjoint of interpolator using precomputed sparse matrix
Yk = full(st.vm * X);		% [K,L]

x = st.K1 * ifft(Yk);	% [K,L] inverse FFT
x = x(1:st.N1,:);	% [N,L] discard zero padding from the end

if size(x,2) > 1
	x = x .* st.sn(:,ones(1,size(x,2)));
else
	x = x .* st.sn;
end
