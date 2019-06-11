 function x = nufft2_adj(X, st)
%function x = nufft2_adj(X, st)
%	SUPERCEDED BY nufft_adj.m
%	Apply adjoint of 2D NUFFT to vector(s) X
%	in:	X	[M,L]
%		st		structure precomputed by nufft2_init()
%	out:	x	[N1,N2,L]
%
%	Copyright 2001-9-17	Jeff Fessler	The University of Michigan

%
%	if no arguments, then run a simple test
%
if nargin < 1
	help(mfilename)

	N1 = 4;
	N2 = 8;
	n_shift = [2.7 3.1];	% random shifts to stress it
	o1 = 2 * pi * [0.0 0.1 0.3 0.4 0.7 0.9]';
	o2 = flipud(o1);
	st = nufft2_init([o1 o2], N1, N2, 4, 6, 2*N1, 2*N2, n_shift, 0, 'best');

	M = length(o1);
	A = zeros(M, N1*N2);
	for n1=1:N1
		for n2=1:N2
			x = zeros(N1,N2); x(n1,n2) = 1;
			jj = n1 + (n2-1) * N1;
			A(:,jj) = nufft2(x, st);
		end
	end

	Aa = zeros(N1*N2, M);
	for ii=1:M
		y = zeros(M,1); y(ii,1) = 1;
		Aa(:,ii) = col(nufft2_adj(y, st));
	end

	printf('max %% = %g', max_percent_diff(A,Aa'))
	clear x
return
end

%	extract attributes from structure
K1 = st.K1;
K2 = st.K2;
if size(X,1) ~= st.M, error size, end

%
%	adjoint of interpolator using precomputed sparse matrix
%
Xk = full(st.p' * X);			% [K1*K2,L]
Xk = reshape(Xk, [K1 K2 size(X,2)]);	% [K1,K2,L]

x = K1 * K2 * ifft2(Xk);		% [K1,K2,L]

% eliminate zero padding from ends
x = x(1:st.N1,1:st.N2,:);		% [N1,N2,L]

if ndims(x) == 2
	x = x .* st.sn;
else
	error '3d not done'
end
