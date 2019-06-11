%{
From smalone@engin.umich.edu Wed Jul 15 23:51:01 1998
To: "Prof. Fessler" <fessler@eecs.umich.edu>
Subject: MATLAB ART Algorithm

Prof. Fessler,

I have two ART algorithms for you to look at.

Like the two PWLS+SOR algorithms, one MATLAB function has G'*W*G as a
precomputed input argument and the other is without that option.

Unlike PWLS+SOR, ART seems to require the full matrix of G'*W*G, not just
the diag. So, in order to precompute, one needs to burn O(n^3) operations.
At least, that's how it seems now.

I still don't understand the ART algorithm in the paper by the Taiwanese
guys.  Maybe we could go over that on Friday.

Brendhan

PS.  Here are my MATLAB functions for ART.

-----ART with GTWG precomputed-----
%}

function [x, r] = art(x, r, G, W, R, y, beta, GTWG)

% -- ART Iterative Routine --
%	Copyright Brendhan Givens, University of Michigan, July 15, 1998
%
%	Details of the basic ART algorithm are described in
%	"Image Reconstruction from Projections", A. C. Kak,
%	4. II. E., Digital Image Processing Techniques,
%	M. Ekstrom (Ed.), Academic Press, Inc., 1984.
%
%	[x, r] = art(x, r, G, W, R, y, beta, GTWG)
%
%	art.m computes one iteration of the overdetermined symmetric
%	positive definite linear system (G'*W*G + beta*R)*x = G'*W*y
%	using the algebraic reconstruction technique with a non-
%	negativity constraint.  Input matrices are assumed to be in
%	sparse format.  Input residual vector must be initialized as
%	r = y - G*x.  GWG is precomputed for efficient computation.
%	x and r are updated.
%
%	input	x	REAL approximation vector
%		r	REAL residual vector
%		G	REAL FULL RANK matrix (sparse)
%		W	REAL diagonal weight matrix (sparse)
%		R	REAL SYMMETRIC penalty matrix (sparse)
%		y	REAL projection data
%		beta	REAL scalar
%		GTWG	REAL SYMMETRIC precomputed matrix where
%			GTWG(i,j) = G(:,i)'*W*G(:,j)
%
%	output	x	REAL NON-NEGATIVE approximation vector
%		r	REAL residual vector

[m, n] = size(G);
w = diag(W);

for j=1:n
	d = ((G(:,j).*w)'*r-beta*R(:,j)'*x)/(GTWG(:,j)'*GTWG(:,j)+2*beta*GTWG(:,j)'*R(:,j)+beta^2*R(:,j)'*R(:,j));
	x_old = x(j);
	x(j) = max(0, x(j) + d);
	r = r + (x_old - x(j))*G(:,j);
end

%% -----ART without precomputed GTWG-----

function [x, r] = art(x, r, G, W, R, y, beta)

%	--	ART Iterative Routine --
%	Brendhan Givens, University of Michigan, July 15, 1998
%
%	Details of the basic ART algorithm are described in
%	"Image Reconstruction from Projections", A. C. Kak,
%	4. II. E., Digital Image Processing Techniques,
%	M. Ekstrom (Ed.), Academic Press, Inc., 1984.
%
%	[x, r] = art(x, r, G, W, R, y, beta)
%
%	art.m computes one iteration of the overdetermined symmetric
%	positive definite linear system (G'*W*G + beta*R)*x = G'*W*y
%	using the algebraic reconstruction technique with a non-
%	negativity constraint.  Input matrices are assumed to be in
%	sparse format.  Input residual vector must be initialized as
%	r = y - G*x.  x and r are updated.
%
%	input	x	REAL approximation vector
%		r	REAL residual vector
%		G	REAL FULL RANK matrix (sparse)
%		W	REAL diagonal weight matrix (sparse)
%		R	REAL SYMMETRIC penalty matrix (sparse)
%		y	REAL projection data
%		beta	REAL scalar
%
%	output	x	REAL NON-NEGATIVE approximation vector
%		r	REAL residual vector

[m, n] = size(G);
w = diag(W);

for j=1:n
	p = (G(:,j).*w);
	pg = G'*p;
	d = (p'*r-beta*R(:,j)'*x)/(pg'*pg+2*beta*pg'*R(:,j)+beta^2*R(:,j)'*R(:,j));
	x_old = x(j);
	x(j) = max(0, x(j) + d);
	r = r + (x_old - x(j))*G(:,j);
end
