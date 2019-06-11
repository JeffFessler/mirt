% test_nufft_sym.m
% test to verify that a signal with conjugate symmetry will produce a real
% spectrum, even with NUFFT approximations.  it does!
% this is useful for the NUFFT-based forward projector.

N = 33;
J = 4;
K = 2*(N-1);
omega = linspace(0, 2*pi, 501)';	% crude spiral:
shift = (N-1)/2;
args = {omega, N, J, K, shift};
G = Gnufft(args);

rng(0)
x = rand((N-1)/2,1) + 1i * rand((N-1)/2,1);
x = [x; 3; conj(flipud(x))];

y = G*x;

ye = dtft(x, omega, shift);

clf, subplot(211)
plot(omega, real(ye), '-', omega, real(y), '--')
subplot(212)
plot(omega, imag(ye), '-', omega, imag(y), '--')
