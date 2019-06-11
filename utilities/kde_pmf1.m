  function [pmf xs] = kde_pmf1(x, varargin)
%|function [pmf xs] = kde_pmf1(x, varargin)
%|
%| Compute an empirical probability mass function (PMF) using interpolation
%| that is similar to a kernel density estimate.
%| This routine is useful for making PMFs (interpolated / smoothed) histograms
%| for computing entropy and other similarity metrics for image registration.
%|
%| in
%|	x	[N 1]	iid realizations of a random variable
%|
%| option
%|	'dx'		bin spacing (default: 1)
%|			use 'silverman' to choose dx based on rule-of-thumb
%|			using inter-quartile range (p48 of 1986 book)
%|	'chat'	0|1	print out intermediate values? (default: 0)
%|
%| out
%|	pmf	[K 1]	nonnegative and sum to 1
%|	xs	[K 1]	sample locations, differing by dx
%|
%| Code assumes kernel support is (-width/2, width/2).
%| Default is a quadratic spline supported on (-3/2,3/2).
%|
%| Copyright 2010-7-4, Jeff Fessler, University of Michigan

arg.dx = 1;
arg.kernel = @kde_pmf1_bspline2;
arg.width = 3; % kernel width, must match kernel!
arg.loop = 'pmf';
arg.nhist = 100;
arg.chat = 0;
arg = vararg_pair(arg, varargin);

if nargin < 1, ir_usage, end

if streq(x, 'test'), kde_pmf1_test, return, end

x = x(:);
[pmf xs] = kde_pmf1_do(x, arg.dx, arg.kernel, arg.width, arg.loop, arg.chat);

if ~nargout && im % plot it if no output requested
	clf, subplot(211)
	[his xh] = hist(x, arg.nhist);
	ax = [min(min(xh), min(xs)) max(max(xh), max(xs))];
	bar(xh, his/numel(x)), title 'histogram'
	axisx(ax)

	subplot(212)
	bar(xs, pmf), title 'KDE PMF'
	axisx(ax)
	colormap('default')
	clear pmf xs
end


% kde_pmf1_do()
function [pmf xs] = kde_pmf1_do(x, dx, kernel, width, loop, chat)

if kernel(-width/2) ~= 0 || kernel(width/2) ~= 0
	fail 'kernel must be zero at +/- width/2 (and beyond)'
end

N = numel(x);

% eqn (3.31) on p 48 of Silverman 1986, which is for gaussian.
% here converted to match the FWHM for quadratic bspline
if streq(dx, 'silverman')
	dx = kde_pmf_width(x, 'chat', chat);
end

kmin = floor(min(x)/dx - width/2) + 1;
kmax = ceil(max(x)/dx + width/2) - 1; % kmin <= k <= kmax
K = kmax - kmin + 1;

switch loop

case 'none' % no loop, but O(NK)
	kk = kmin:kmax;
	pmf = sum(kernel(outer_sum(x/dx, -kk)), 1)';

case 'data' % loop over data N
	pmf = zeros(K,1);
	for nn=1:N
		xn = x(nn);
		k1 = floor(xn/dx - width/2) + 1;
		k2 = ceil(xn/dx + width/2) - 1;
		kk = (k1:k2)';
		pmf(kk-kmin+1) = pmf(kk-kmin+1) + kernel(xn/dx - kk);
	end

case 'pmf' % loop over PMF bins K
	pmf = zeros(K,1);
	for ik=1:K
		kk = kmin + ik - 1;
		xmin = ((kk - width/2)) * dx;
		xmax = ((kk + width/2)) * dx;
		good = xmin < x & x < xmax;
		pmf(ik) = sum(kernel(x(good)/dx - kk));
	end

case 'kernel' % loop over kernel width
	fail 'does not work because multiple identical k values per'
	pmf = zeros(K+0,1);
	k1 = floor(x/dx - width/2) + 1;
	for ll=0:ceil(width)-1
		kk = k1 + ll;
		pmf(kk-kmin+1) = pmf(kk-kmin+1) + kernel(x/dx - kk);
	end
%	if any(pmf(K+1:end)), keyboard, end
%	pmf = pmf(1:K);

otherwise
	fail('unknown loop type "%s"', loop)
end
xs = (kmin:kmax)' * dx;
pmf = pmf / N;


% quadratic B-spline
% kde_pmf1_bspline2()
function y = kde_pmf1_bspline2(x)
y = (3/4 - x.^2) .* (abs(x) < 1/2) + ...
	(abs(x) - 3/2).^2 / 2 .* (1/2 <= abs(x) & abs(x) < 1.5);


% test routine
function kde_pmf1_test

rng(0)
n = 300;
sig = 5;
x = sig * randn(n, 1);

if 1 % check bspline2
	t = linspace(-2,2,101);
	dt = t(2) - t(1);
	b2 = kde_pmf1_bspline2(t);
	gauss = @(x,s) 1/sqrt(2*pi*s^2) * exp(-x.^2 / 2 / s^2);
	if im
		clf, subplot(211)
		fwhm = 3 - sqrt(3); % FWHM of bspline2, about 1.28
		s = fwhm / sqrt(log(256)); % sig of gauss with matching FWHM
		gauss0 = @(x) gauss(x,s)/gauss(0,s)*3/4; % scaled
		plot(t, b2, '-', t, gauss0(t), '--', ...
			fwhm/2, gauss0(fwhm/2), 'o')
		legend('b2', 'gauss', 'FWHM/2')
		subplot(212)
		plot(t, diffc(b2) / dt)
	prompt
	end
end

dx = 1.8;
[pmf xs] = kde_pmf1(x, 'dx', dx);

[pmf xs] = kde_pmf1(x, 'dx', 'silverman', 'chat', 1);
if 1 % check other loop methods
	loop_list = {'none', 'pmf', 'data'}; % 'kernel'
	for ii=1:numel(loop_list)
		[pmfn xsn] = kde_pmf1(x, 'dx', 'silverman', ...
			'loop', loop_list{ii});
		equivs(pmf, pmfn)
		jf_equal(xs, xsn)
	end
end

equivs(sum(pmf), 1)
h = (1 + log(2*pi*sig^2)) / 2; % differential entropy for gaussian
pr h
Hf = @(p) sum(-p(p>0) .* log(p(p>0)));
H = Hf(pmf);
pr H
pr H + dx * log(dx)

if im
	im clf, pl = @(t) subplot(130 + t);
	pl(1)
	[his xh] = hist(x, 20);
	plot(xh, his/n, '.-')

	pl(2)
	plot(xs, pmf, '.-')

	dx = linspace(0.2, 2.8, 31);
	nd = numel(dx);
	for id=1:nd
		pmf = kde_pmf1(x, 'dx', dx(id));
		Hv(id) = Hf(pmf); % H for several dx values
	end

	pl(3)
	plot(dx([1 end]), [h h], '-', ...
		dx, Hv, 'o', ...
		dx, Hv + log(dx), '+') % H + log(dx) approximately is h
	legend('h', 'H', 'H + log(dx)')
end
