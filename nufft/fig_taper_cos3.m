% fig_taper_cos3
% compare 3 different interpolation kernels:
%	min-max kernel with uniform scaling factors
%	min-max kernel with optimized scaling factors
%	cos^3-tapered dirichlet kernel of magnusson

J = 6;
N = 1000;
K = 2 * N;
k = linspace(-J/2,J/2,101)';
fm = diric(2*pi*k/J, J) .* cos((2*pi*k/J)/2).^3;	% tapered
fn = nufft1_kernel(k, N, J, K, 1, 0);			% min-max unif
fo = nufft1_kernel(k, N, J, K, 'best');
fo = fo / max(fo);

%fm = 0*fm;
%fn = nufft1_kernel(k, N, J, 9*N, 1, 0);
%fo = nufft1_kernel(k, N, J, 2*N, 1, 0);

plot(k, fn, '-', k, fm, '--', k, fo, '-.')
axisy(-0.3,1.1)
legend('min-max', 'Magnusson', 'optimized')
title(sprintf('J=%d K/N=%g', J, K/N))
