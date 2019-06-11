From smalone@engin.umich.edu Wed Jul 15 13:37:52 1998
To: "Jeffrey A. Fessler" <fessler@eecs.umich.edu>
Subject: Re: pwls_sor

Prof. Fessler,

I'm sending along two of my pwls_sor routines.  Both compute one iteration
of the routine.  One uses precomputed diag[G'*W*G].  The other does not.
I could combine the two, if you want, using the nargin command to
determine which one to use.  But I'm assuming your mostly interested in
the routine that uses the precomputed vector.  Let me know if you need any
changes to the code.

Brendhan

----pwls_sor with precomputed diag[G'*W*G]---

function [x, r] = pwls_sor(x, r, G, W, R, y, beta, GTWG)

%    --	PWLS+SOR Iterative Routine --
%	Copyright Brendhan Givens, University of Michigan, July 9, 1998
%
%	Details of this algorithm are described in
%	"Penalized Weighted Least-Squares Image Reconstruction
%	for Positron Emission Tomography", Jeffrey A. Fessler,
%	IEEE Transactions on Medical Imaging, Vol. 13, No. 2,
%	June 1994.
%
%	[x, r] = pwls_sor(x, r, G, W, R, y, beta, GTWG)
%
%	pwls_sor.m computes one iteration of the overdetermined symmetric
%	positive definite linear system (G'*W*G + beta*R)x = G'*W*y using
%	a sequential minimization procedure with a non-negativity
%	constraint.  Input matrices are assumed to be in sparse format.
%	Input residual vector must be initialized as r = y - G*x.  The
%	diag[G'*W*G] is precomputed for faster iterations.
%
%	input	x	REAL approximation vector
%		r	REAL residual vector
%		G	REAL FULL RANK matrix (sparse)
%		W	REAL diagonal weight matrix (sparse)
%		R	REAL SYMMETRIC penalty matrix (sparse)
%		y	REAL projection data
%		beta	REAL scalar
%		GTWG	REAL SYMMETRIC POSITIVE INDEFINITE vector
%			where elements gtwg(j) = G(:,j)'*W*G(:,j)
%
%	output	x	REAL NON-NEGATIVE approximation vector
%		r	REAL residual vector

[m,n] = size(G);			% initialization
w = diag(W);

for j = 1:n				% begin iteration
   d = ((G(:,j).*w)'*r - beta*R(:,j)'*x)/(GTWG(j) + beta*R(j,j));
   x_old = x(j);
   x(j) = max(0, x(j) + d);		% update with non-negativity
constraint
   r = r + (x_old - x(j))*G(:,j);	% update residual
end

----pwls_sor without precompute----

function [x, r] = pwls_sor(x, r, G, W, R, y, beta)

%    --	PWLS+SOR Iterative Routine --
%	Brendhan Givens, University of Michigan, July 9, 1998
%
%	Details of this algorithm are described in
%	"Penalized Weighted Least-Squares Image Reconstruction
%	for Positron Emission Tomography", Jeffrey A. Fessler,
%	IEEE Transactions on Medical Imaging, Vol. 13, No. 2,
%	June 1994.
%
%	[x, r] = pwls_sor(x, r, G, W, R, y, beta)
%
%	pwls_sor.m computes one iteration of the overdetermined symmetric
%	positive definite linear system (G'*W*G + beta*R)x = G'*W*y using
%	a sequential minimization procedure with a non-negativity
%	constraint.  Input matrices are assumed to be in sparse format.
%	Input residual vector must be initialized as r = y - G*x.
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

[m,n] = size(G);			% initialization
w = diag(W);

for j = 1:n				% begin iteration
   p = (G(:,j).*w);
   d = (p'*r - beta*R(:,j)'*x)/(p'*G(:,j) + beta*R(j,j));
   x_old = x(j);
   x(j) = max(0, x(j) + d);		% update with non-negativity
constraint
   r = r + (x_old - x(j))*G(:,j);	% update residual
end



On Wed, 15 Jul 1998, Jeffrey A. Fessler wrote:

>
> please email me your pwls_sor.m routine - the one that is just
> the update without the loop wrapper liked we discussed.  there is
> a nuc. eng. student who wants to use it.
> in fact, i think i'll put a collection of such routines on my web page!
> (with your name on this one of course).
> jf
>


