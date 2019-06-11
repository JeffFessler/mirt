 function [err,f1,f2] = nufft2_err_mm(om, N1, N2, J1, J2, K1, K2, alpha, beta)
%function [err,f1,f2] = nufft2_err_mm(om, N1, N2, J1, J2, K1, K2, alpha, beta)
% Compute worst-case error for each input frequency for min-max 2D NUFFT.
% in:
%	om	[M,2]	digital frequency omega in radians
%	N1,N2		signal length
%	J1,J2		# of neighbors used per frequency location
%	K1,K2		FFT size (should be > N1)
%	alpha	[L,1]	Fourier series coefficients of scaling factors
%	beta		scale gamma=2pi/K by this in Fourier series
%				typically is K/N (me) or 0.5 (Liu)
% out:
%	err	[M,1]	worst-case error over unit-norm signals
%
% Copyright 2001-12-7, Jeff Fessler, The University of Michigan

% if no arguments, give an example
if nargin < 4
	help(mfilename)
	N1 = 1; K1 = 2*N1; J1 = 7; gam1 = 2*pi/K1;
	N2 = 1; K2 = 2*N2; J2 = 6; gam2 = 2*pi/K2;
	alpha = 1;
	alpha = 'best';
	[err,f1,f2] = nufft2_err_mm('all', N1, N2, J1, J2, K1, K2, alpha);
	mesh(f1, f2, err)
	plot(f1, err)
	xlabel '\omega_1 / \gamma', ylabel 'E_{max}(\omega_1,\omega_2)'
	clear err
return
end

if ~isvar('alpha') || isempty(alpha)
	alpha = [1];	% default Fourier series coefficients of scaling factors
end
if ~isvar('beta') || isempty(beta)
	beta = 0.5;	% default is Liu version for now
end
alpha1 = alpha;
alpha2 = alpha;
beta1 = beta;
beta2 = beta;

%
% trick to look at all relevant om's
%
if ischar(om)
	if ~streq(om, 'all'), error 'unknown om argument', end
	gam1 = 2*pi/K1;
	gam2 = 2*pi/K2;
	[f1,f2] = ndgrid(linspace(0,1,31), linspace(0,1,33));
	om = [gam1*f1(:) gam2*f2(:)];
end

%
% see if 'best' alpha is desired
%
if ischar(alpha)
	if ~streq(alpha, 'best'), error 'unknown alpha argument', end
	[alpha1, beta1, ok] = nufft_best_alpha(J1, 0, K1/N1);
	if ~ok, error 'unknown J,K/N', end
	[alpha2, beta2, ok] = nufft_best_alpha(J2, 0, K2/N2);
	if ~ok, error 'unknown J,K/N', end
end

tol = 0;
T1 = nufft_T(N1, J1, K1, tol, alpha1, beta1);		% [J,J]
T2 = nufft_T(N2, J2, K2, tol, alpha2, beta2);		% [J,J]
r1 = nufft_r(om(:,1), N1, J1, K1, alpha1, beta1);	% [J,M]
r2 = nufft_r(om(:,2), N2, J2, K2, alpha2, beta2);	% [J,M]
T = kron(T2, T1);	% J1*J2 x J1*J2

% kron for each om
M = size(om,1);
r = zeros(J1*J2,M);
for ii=1:M
	r(:,ii) = kron(r2(:,ii), r1(:,ii));
end

%
% worst-case error at each frequency
%
Tr = T * r;			% [J,M]
err = sum(conj(r) .* Tr).';	% [M,1]
err = reale(err);
err = min(err, 1);
err = sqrt(1 - err);

if isvar('f1')
	err = reshape(err, size(f1));
end
