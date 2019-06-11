 function centers = highrate_centers(data, L, M)
%function centers = highrate_centers(data, L, M)
%| According to high-rate scalar quantization theory (Gersho, T-IT, 1979),
%| the density of MMSE quantization cell centroids should be f^{k/(k+2)},
%| for k-dimensional data distributed with pdf f.
%| This m-file designs L centroids for scalar data (k=1) using this principle.
%| M is the number of histogram bins used to approximate the data pdf. 
%| out: centers is [L,1]
%| Copyright 2004-7-7, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if streq(data, 'test'), highrate_centers_test, return, end

if ~isvar('M'), M=100; end

data = data(:);
[wt data] = hist(data, M);
dens = wt .^ (1/3); % from high-rate scalar quantization theory
cdf = cumsum(dens / sum(dens));
[cdf ii] = unique(cdf);
m1 =  [1:M]-0.5;
m1 = m1(ii);
uu = ([1:L]-0.5)/L;
m2 = interp1(cdf, m1, uu, 'pchip');
%plot(cdf, m1, '.', uu, m2, 'o')
m2 = max(m2,1);
%try
centers = data(round(m2));
%catch
%minmax(m2)
%keyboard
%end

function highrate_centers_test
rng(0)
N = 10^6;
g = rand(N, 1);
f = 2 ./ (2 - g); fun = @(f) 2./f.^2 .* (1 <= f & f <= 2); % pdf on [1,2]
%f = 3+2*g; fun = @(f) 1/2 .* (3 <= f & f <= 5); % pdf on [3,5]
im plc 2 2
im subplot 1
if 1
	hist(f, 100)
	hold on
	u = linspace(1,2,100); du = u(2) - u(1);
	plot(u, fun(u)*N*du, 'g-')
	hold off
	xlabelf '$f$', ylabelf 'histogram'
	xtick([1 2]), ytick([0 2e4])
	axis([1 2 0 2e4])
end
L = 2^4; % even this is high enough for this example!
l1 = [1:L]';
rl_compand = 2 ./ (2 - (l1-1/2)/L); % compand levels for 2/f^2
%rl_ideal = 3 + (l1-1/2)/L * (5-3); % ideal levels for U(3,5)
rl_high = highrate_centers(f, L, 2^9);

% todo: run lloyd-max to find "ideal" and compare!
rl_ideal = lloyd_max_hist(f, rl_compand, 10^3);

im subplot 2
plot(	l1, rl_high, '-+', ...
	l1, rl_compand, '-o', ...
 	l1, rl_ideal, '-*')
legend('high rate', 'companding', 'lloyd-max', 'location', 'northwest')
axis([1 L 1 2])
xlabelf '$l$', ylabelf '$r_l$'
xtick([1 L]), ytick([1 2])

% ir_savefig cw highrate_centers
