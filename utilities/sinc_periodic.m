  function x = sinc_periodic(t, K)
%|function x = sinc_periodic(t, K)
%| periodic sinc function, obtained by replicates of a sinc:
%| x(t) = \sum_l sinc(t - l K) = ... = sin(pi*t) / tan(pi*t/K) / K for K even.
%| This function is bandlimited and its samples are an impulse train.
%| It is closely related to the Dirichlet function diric() for odd K,
%| but it differs for even K.
%| Copyright 2003-11-2, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if streq(t, 'test'), sinc_periodic_test, return, end

if ~rem(K,2) % even
	d = tan(pi*t/K);
	j = abs(d) > 1e-12;
	x = ones(size(t));
	t = t(j);
	d = d(j);
	x(j) = sin(pi*t) ./ d / K;
else
	x = nufft_diric(t, K, K, 1);
%	xo = diric(2*pi*t/K,K); % would require matlab's signal toolbox
%	max_percent_diff(xo, x)
end


%
% self test
%
function sinc_periodic_test

Klist = [4 5];
im clf, pl=220;
for kk=1:2
	K = Klist(kk);
	n = [0:(4*K)]';
	t = linspace(0,4*K,401)';
	x = @(t,K) sinc_periodic(t, K);
%	y = @(t,K) diric(2*pi*t/K,K); % would require matlab's signal toolbox
	y = @(t,K) nufft_diric(t,K,K,1);
	if im
		subplot(pl+kk+0)
		plot(t, x(t,K), '-', n, x(n,K), 'o')
		titlef('Sinc-Periodic K=%d', K)
		axis([0 4*K -1 1]), xtick([0:4]*K), grid
		subplot(pl+kk+2)
		plot(t, y(t,K), '-', n, y(n,K), 'o')
		titlef('Dirichlet K=%d', K)
		axis([0 4*K -1 1]), xtick([0:4]*K), grid
		ylabelf('K=%d', K)
	end
end
