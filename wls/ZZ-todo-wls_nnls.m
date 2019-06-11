function [x,w] = nn_rls(E,R,f,tol)
%	Non-negative REGULARIZED least-squares.
%	X = NNLS(A,b) returns the vector X that solves A*x = b
%	in a least squares sense, subject to x >= 0.
%
%	A default tolerance of TOL = MAX(SIZE(A)) * NORM(A,1) * EPS
%	is used for deciding when elements of X are less than zero.
%	This can be overridden with X = NNLS(A,b,TOL).
%
%	[X,W] = NNLS(A,b) also returns dual vector W where
%	w(i) < 0 where x(i) = 0 and w(i) = 0 where x(i) > 0.
%	L. Shure 5-8-87
%	Revised, 12-15-88 LS.
%	Copyright (c) 1987-89 by the MathWorks, Inc.
%	modifed by W. Ahmad and J. Fessler

% Reference:
%  Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.

% initialize variables
if nargin < 4 % Changed it from 3 to 4 (no. of input arguments)
	tol = eps*norm(E,1)*max(size(E));
end
[m,n] = size(E);
P = zeros(1,n);
Z = 1:n;
x = P';
ZZ=Z;
w = E'*(f-E*x)-R*x;

n_iter = 0;
% outer loop to put variables into set to hold positive coefficients
while any(Z) & any(w(ZZ) > tol)
	n_iter = n_iter + 1;
	[wt,t] = max(w(ZZ));
	t = ZZ(t);
	P(1,t) = t;
	Z(t) = 0;
	PP = find(P);
	ZZ = find(Z);
	nzz = size(ZZ);
%	EP(1:m,PP) = E(:,PP);
%	EP(:,ZZ) = zeros(m,nzz(2));
%	z = inv(EP'*EP+R)*EP'*f;
%%	z=pinv(EP)*f;
	MM = E' * E + R;	% K
	MP = MM(PP,PP);
	z = zeros(n,1);
	z(PP) = inv(MP) * E(:,PP)' * f;	% K
	z(ZZ) = zeros(nzz(2),nzz(1));

% inner loop to remove elements from the positive set which no longer belong
	while any((z(PP) <= tol))
%		disp(find(z(PP) <= tol))
		QQ = find((z <= tol) & P');
		alpha = min(x(QQ)./(x(QQ) - z(QQ)));
		x = x + alpha*(z - x);
		ij = find(abs(x) < tol & P' ~= 0);
		Z(ij)=ij';
		P(ij)=zeros(1,max(size(ij)));
		PP = find(P);
		ZZ = find(Z);
		nzz = size(ZZ);
%		EP(1:m,PP) = E(:,PP);
%		EP(:,ZZ) = zeros(m,nzz(2));
%		z = inv(EP'*EP+R)*EP'*f;
%%		z = pinv(EP)*f;
		MM = E' * E + R;	% K
		MP = MM(PP,PP);
		z = zeros(n,1);
		z(PP) = inv(MP) * E(:,PP)' * f;	% K
		z(ZZ) = zeros(nzz(2),nzz(1));
	end
        x = z;
	w = E'*(f-E*x)-R*x;
end
%	plot(x)
% n_iter
