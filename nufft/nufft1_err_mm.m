 function [err, sn, T1] = nufft1_err_mm(om, N1, J1, K1, type, alpha, beta)
%function [err, sn, T1] = nufft1_err_mm(om, N1, J1, K1, type, alpha, beta)
% Compute worst-case error for each input frequency for min-max 1D NUFFT.
% in:
%	om	[M,1]	digital frequency omega in radians
%	N1		signal length
%	J1		# of neighbors used per frequency location
%	K1		FFT size (should be > N1)
%	type		'sinc' 'diric' 'qr'
%	alpha	[L,1]	Fourier series coefficients of scaling factors
%			trick: or, "sn" if length N1
%	beta		scale gamma=2pi/K by this in Fourier series
%			typically is K/N (me) or 0.5 (Liu)
% out:
%	err	[M,1]	worst-case error over unit-norm signals
%	sn	[N,1]	scaling factors corresponding to alpha,beta
%	T1	[J,J]	T matrix
%
% Copyright 2001-12-7, Jeff Fessler, The University of Michigan

% if no arguments, give an example
if nargin < 4
	help(mfilename)
	N = 100; K = 2*N; gam = 2*pi/K;
	J = 14;
	om = gam * linspace(0,1,101);
	[alpha, beta, ok] = nufft_best_alpha(J, 2, K/N);
	if ~ok, alpha = []; beta = 0.5; end
	errd = nufft1_err_mm(om, N, J, K, 'diric');
	errs = nufft1_err_mm(om, N, J, K, 'sinc');
	errq = nufft1_err_mm(om, N, J, K, 'qr');
	semilogy(om/gam, errs, 'g-x', om/gam, errd, 'y-+', om/gam, errq, 'c-o')
	xlabel '\omega / \gamma', ylabel 'E_{max}(\omega)'
	legend('Tr sinc', 'Tr diric', 'QR approach')
return
end

if ~isvar('type') || isempty(type),	type = 'sinc'; end
if ~isvar('alpha') || isempty(alpha)
	alpha = [1];	% default Fourier series coefficients of scaling factors
end
if ~isvar('beta') || isempty(beta)
	beta = 0.5;	% default is Liu version for now
end

use_qr = false;
if streq(type, 'sinc')
	use_true_diric = false;
elseif streq(type, 'diric')
	use_true_diric = true;
elseif streq(type, 'qr')
	use_qr = true;
else
	error 'unknown type'
end


%
% see if 'best' alpha is desired
%
if ischar(alpha)
	if streq(alpha, 'uniform')
		alpha = [1];
		beta = 0.5;
	else
		if streq(alpha, 'best')
			L = 0;
		elseif streq(alpha, 'best,L=1')
			L = 1;
		elseif streq(alpha, 'best,L=2')
			L = 2;
		else
			error 'unknown alpha argument'
		end
		[alpha, beta, ok] = nufft_best_alpha(J1, L, K1/N1);
		if ~ok
			tmp = 'optimal alpha unknown for J=%d, K/N=%g, L=%d';
			warning(sprintf(tmp, J1, K1/N1, L))
			sn = ones(N1,1);
			err = nan;
			return
		end
	end
end

%
% if requested, return corresponding scaling factors too
%
if length(alpha) == N1	% trick: special way to give explicit sn's
	sn = alpha;
	if ~use_qr, error 'give sn only allowed for QR version', end
elseif nargout > 1 || use_qr
	sn = nufft_scale(N1, K1, alpha, beta);
end

%
% QR approach to error
%
if use_qr
	n = [0:N1-1]' - (N1-1)/2;
	[nn, jj] = ndgrid(n, 1:J1);
	gam1 = 2*pi/K1;
	C = exp(1i * gam1 * nn .* jj) / sqrt(N1);
	S = spdiag(sn, 'nowarn');
	A = S' * C;
	[Q R] = qr(S' * C, 0);	% [N,J] compact QR decomposition

	do = col(om - gam1*nufft_offset(om, J1, K1));
	Db = exp(1i * n * do') / sqrt(N1);	% [N,M]
	err = Db - Q * (Q' * Db);
	err = sqrt(sum(abs(err).^2,1))';	% [M]
return
end

tol = 0;
T1 = nufft_T(N1, J1, K1, tol, alpha, beta, use_true_diric);	% [J,J]
r1 = nufft_r(om, N1, J1, K1, alpha, beta, use_true_diric);	% [J,M]

%
% worst-case error at each frequency
%
Tr1 = T1 * r1;			% [J,M]
err = sum(conj(r1) .* Tr1).';	% [M,1]
err = min(real(err), 1);
err = sqrt(1 - err);		% caution: this "1 -" may cause numerical error
