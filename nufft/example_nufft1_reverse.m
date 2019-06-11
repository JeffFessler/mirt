% "example_nufft1_reverse.m"
% This m-file is an example of applying the NUFFT method "in reverse,"
% meaning that you have a nonuniformly-sampled signal spectrum in hand
% and you want to try to "compute its inverse Fourier transform"
% to get uniform set of signal samples.  In the nomenclature of the 1999
% SIAM J Sci. Comput. paper by Nguyen and Liu, this is "Problem 1."
% (It is somewhat related to MRI reconstruction by gridding.)
% Before doing this, you should ask yourself if that is *really*
% what you want to do.  In most inverse problems, I think it is NOT
% what we should do, as argued in the T-SP paper on this method.
% Instead, we should solve the inverse problem *iteratively*, using
% the NUFFT method (and its adjoint) each iteration.
% Furthermore, the min-max optimality of the NUFFT method was only
% established for the "forward direction:" going from uniform signal
% samples to nonuniform frequency samples.
% (You can replace time and frequency in the above discussion.)
% Anyway, let's try it here and see how it works...

%
% synthesize some spectral data that is nonuniformly spaced
%
No = 2^7;				% number of frequency samples
rng(0)
om = 2*pi*sort(rand(No,1)-0.5);		% random frequency samples on [-pi,pi ]
%om = 2*pi*[-No/2:No/2-1]'/No;		% test with uniform samples

% a spectrum with periodic components, hopefully visible in other domain
X = inline('2 + 2*sin(om*10) + 4*cos(15*om) + 4*cos(20*om)', 'om');
Xm = X(om);

%
% go to the time domain the "exact" (slow) DTFT way.
% this implements essentially equation (1) in Nguyen and Liu (1999)
%
N1 = 2^6;		% # of time samples
n = [0:N1-1]'-N1/2;	% time sample locations
xd = dtft2_adj(Xm, [om 0*om], N1, 1, [N1/2 0]);

%
% now do it the fast NUFFT way
%
if 1 || ~isvar('st')	% create NUFFT structure
	J = 5;		% interpolation neighborhood
	K1 = N1*2;	% two-times oversampling
	st = nufft_init(om, N1, J, K1, N1/2, 'minmax:kb');
end

xn = nufft_adj(Xm, st);	% call ADJOINT to go "in reverse"

%
% compare slow exact to fast NUFFT
%
figure(1), clf, pl=220;
subplot(pl+1)
oo = [-200:200]'/400*2*pi;
plot(om, Xm, '.', oo, X(oo), '-')
xlabel \omega, ylabel X(\omega), title 'Spectrum and samples'

subplot(pl+2)
plot(n, real(xd), 'c.-', n, imag(xd), 'y.-')
xlabel n, ylabel x[n], title 'Exact "nonuniform FT"'

subplot(pl+3)
plot(n, real(xn), 'c.-', n, imag(xn), 'y.-')
xlabel n, ylabel x[n], title 'Fast NUFFT adjoint'

subplot(pl+4)
plot(n, real(xn-xd), 'c.-', n, imag(xn-xd), 'y.-')
xlabel n, ylabel x[n], title 'Approximation Error (very small!)'


%
% ok fine, so if you look at figure 1 you see that
% the NUFFT gives a great approximation to the "exact" formula
% but is the result *really* what you wanted?
% It seems that ideally the result should be 4 spikes with
% no background junk.  If we solved "Problem 5" in Nguyen and Liu,
% then we should get that!  This is how we have approached the
% MRI reconstruction problem, as described in ../mri.
%
% For iterative 2D version, see ../example/mri_example.m

return

%
% here is a different strategy: interpolate the nonuniform data
% onto a uniform grid, then simply take the inverse FFT.
% this didn't work so well because it is a lousy gridding method.
% todo: need to put a better gridding method here!
%
ok = [-N1/2:(N1/2-1)]'/N1*2*pi;
Xg = interp1(om, Xm, ok, 'linear', 'extrap');
xg = fftshift(ifft(fftshift(Xg)));

%xg = xg ./ nufft_sinc(n/N1).^2; % post-compensate?

figure(2), clf, pl=220;
subplot(pl+1)
plot(n, real(xg), 'c.-', n, imag(xg), 'y.-')
xlabel n, ylabel x[n], title 'Gridding "reconstruction"'

Xg = dtft1(xg, om, N1/2);
Xd = dtft1(xd/No, om, N1/2);

subplot(pl+2)
plot(om, real(Xg), '.', om, real(Xm), '-')
xlabel \omega, ylabel X(\omega), title ''

subplot(pl+3)
plot(om, real(Xd), '.', om, real(Xm), '-')
xlabel \omega, ylabel X(\omega), title ''
