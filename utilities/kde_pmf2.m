 function [pmf xs ys] = kde_pmf2(x, y, varargin)
%function [pmf xs ys] = kde_pmf2(x, y, varargin)
%|
%| Compute a 2D empirical joint probability mass function (PMF)
%| using interpolation that is similar to a kernel density estimate.
%| This joint PMF (interpolated / smoothed histogram) is useful for
%| computing joint entropy and other similarity metrics for image registration.
%|
%| in
%|	x	[N 1]
%|	y	[N 1]	iid realizations of random variable pairs (x,y)
%|
%| option
%|	'dx'		bin spacing (default: 1)
%|			use 'silverman' to choose dx based on rule-of-thumb
%|			using inter-quartile range (p48 of 1986 book)
%|	'dy'		bin spacing (default: dx)
%|	'loop'		method; default is to try 'mex' if available, else 'pmf'
%|	'chat'	0|1	print out intermediate values? (default: 0)
%|
%| out
%|	pmf	[Kx Ky]	nonnegative and sum to 1
%|	xs	[Kx 1]	sample locations, differing by dx
%|	ys	[Ky 1]	sample locations, differing by dy
%|
%| Code assumes kernel support is (-width/2, width/2).
%| Default is a quadratic spline supported on (-3/2,3/2).
%|
%| Copyright 2010-07-31, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), kde_pmf2_test, return, end
if nargin < 2, ir_usage, end

arg.dx = 1;
arg.dy = []; % dx below
arg.kernel = @kde_pmf2_bspline2;
arg.widthx = 3; % kernel width, must match kernel!
arg.widthy = []; % widthx below
arg.loop = ''; % see below
arg.chat = 0;
arg = vararg_pair(arg, varargin);

if isempty(arg.loop)
	if has_mex_jf
		arg.loop = 'mex';
	else
		arg.loop = 'pmf';
	end
end
if isempty(arg.dy)
	arg.dy = arg.dx;
end
if isempty(arg.widthy)
	arg.widthy = arg.widthx;
end

x = x(:);
y = y(:);

% eqn (3.31) on p 48 of Silverman 1986, which is for gaussian.
% here converted to match the FWHM for quadratic bspline
if streq(arg.dx, 'silverman')
	arg.dx = kde_pmf_width(x, 'chat', arg.chat);
end
if streq(arg.dy, 'silverman')
	arg.dy = kde_pmf_width(y, 'chat', arg.chat);
end

if streq(arg.loop, 'mex')
	if ~isequal(arg.kernel, @kde_pmf2_bspline2) ...
		|| arg.widthx ~= 3 || arg.widthy ~= 3
		fail 'mex only supports bspline2 with width=3'
	end
	[pmf xs ys] = jf_mex('kde,b2', single([x y]), single([arg.dx arg.dy]));

else
	[pmf xs ys] = kde_pmf2_do(x, y, ...
		arg.dx, arg.dy, arg.widthx, arg.widthy, ...
		arg.kernel, arg.kernel, arg.loop, arg.chat);
end

if ~nargout && im
	im plc 1 2
	[his cen] = jf_histn([x(:) y(:)]);
	im(1, cen{1}, cen{2}, his, 'histogram')
	im(2, xs, ys, pmf, 'pmf')
	clear pmf
end


% kde_pmf2_do()
function [pmf xs ys] = kde_pmf2_do(x, y, dx, dy, wx, wy, kernx, kerny, ...
	loop, chat)

if kernx(-wx/2) ~= 0 || kernx(wx/2) ~= 0
	fail 'kernel must be zero at +/- width/2 (and beyond)'
end
if kerny(-wy/2) ~= 0 || kerny(wy/2) ~= 0
	fail 'kernel must be zero at +/- width/2 (and beyond)'
end

N = numel(x);
if N ~= numel(y), fail('x and y must be same size'), end

kmin = floor(min(x)/dx - wx/2) + 1;
kmax = ceil(max(x)/dx + wx/2) - 1; % kmin <= k <= kmax
Kx = kmax - kmin + 1;
xs = (kmin:kmax)' * dx;
x = x - dx * (kmin-1); % wlog, shift so that kmin=1 (for matlab)
%if chat, pr [kmin kmax], end

kmin = floor(min(y)/dy - wy/2) + 1;
kmax = ceil(max(y)/dy + wy/2) - 1; % kmin <= k <= kmax
Ky = kmax - kmin + 1;
ys = (kmin:kmax)' * dy;
y = y - dy * (kmin-1); % wlog, shift so that kmin=1 (for matlab)
%if chat, pr [kmin kmax], end
clear kmin kmax

switch loop

case 'none' % no loop, but O(NK) and huge memory use!
	kx = 1:Kx;
	ky = 1:Ky;
	tmp1 = kernx(outer_sum(x/dx, -kx)); % [N Kx]
	tmp2 = kerny(outer_sum(y/dy, -ky)); % [N Ky]
	pmf = tmp1' * tmp2; % [Kx Ky]

case 'data' % loop over data N
	pmf = zeros(Kx,Ky);
	for nn=1:N
		xn = x(nn);
		yn = y(nn);
		k1 = floor(xn/dx - wx/2) + 1;
		k2 = ceil(xn/dx + wx/2) - 1;
		kx = (k1:k2)';
		k1 = floor(yn/dy - wy/2) + 1;
		k2 = ceil(yn/dy + wy/2) - 1;
		ky = (k1:k2);
		pmf(kx,ky) = pmf(kx,ky) + ...
			kernx(xn/dx - kx) * kerny(yn/dy - ky); % outer prod
	end

case 'pmf' % loop over PMF bins K
	pmf = zeros(Kx,1);
	for ky=1:Ky
	 for kx=1:Kx
		xmin = ((kx - wx/2)) * dx;
		xmax = ((kx + wx/2)) * dx;
		good = xmin < x & x < xmax;
		ymin = ((ky - wy/2)) * dy;
		ymax = ((ky + wy/2)) * dy;
		good = (ymin < y & y < ymax) & good;
		tmp = kernx(x(good)/dx - kx) .* kerny(y(good)/dy - ky);
		pmf(kx,ky) = sum(tmp);
	 end
	end

otherwise
	fail('unknown loop type "%s"', loop)
end
pmf = pmf / N;


% quadratic B-spline
% kde_pmf2_bspline2()
function y = kde_pmf2_bspline2(x)
y = (3/4 - x.^2) .* (abs(x) < 1/2) + ...
	(abs(x) - 3/2).^2 / 2 .* (1/2 <= abs(x) & abs(x) < 1.5);


% test routine
function kde_pmf2_test

rng(0)
n = 500;
sig = 5;
rho = 0.9; % correlated gaussian
Cov = sig^2 * [1 rho; rho 1];

%h1 = (1 + log(2*pi*sig^2)) / 2; % differential entropy for 1D gaussian
h = (size(Cov,1) + log(det(2*pi*Cov))) / 2; % differential entropy for gaussian
pr h

tmp = randn(n, 2) * sqrtm(Cov);
x = tmp(:,1);
y = tmp(:,2);

if 0 % uniform
	x = 10 * (rand(n,1) - 0.5);
	y = 10 * (rand(n,1) - 0.5);
end

if 1 % check h
	dx = 2.4;
	dy = 2.8;
	[pmf xs ys] = kde_pmf2(x, y, 'dx', dx, 'dy', dy);
	equivs(sum(pmf(:)), 1)
	Hf = @(p) sum(-p(p>0) .* log(p(p>0)));
	H = Hf(pmf);
	pr H
	pr H + log(dx*dy)
end

args = {x, y, 'dx', 'silverman'};
[pmf xs ys] = kde_pmf2(args{:}, 'chat', 1);

if 1 % check other loop methods
	loop_list = {'none', 'pmf', 'data'};
	if has_mex_jf
		loop_list{end+1} = 'mex';
	end
	for ii=1:numel(loop_list)
		cpu etic
		[pmfn xsn ysn] = kde_pmf2(args{:}, 'loop', loop_list{ii});
		equivs(pmf, pmfn, 'thresh', 2e-6)
		equivs(xs, xsn)
		equivs(ys, ysn)
		cpu('etoc', sprintf('%s: %4s', mfilename, loop_list{ii}))
	end
end

if 0 && has_mex_jf
	[pmfm xsm ysm] = kde_pmf2(args{:}, 'loop', 'mex');
%	im(xsm, ysm, pmfm)
	equivs(xsm, xs)
	equivs(ysm, ys)
	equivs(pmf, pmfm)
end

if im
	im plc 1 3
	[hist cent] = jf_histn([x y], 'nbin', 30);
	im(1, cent{1}, cent{2}, hist/n)
	im(2, xs, ys, pmf)

	dx = linspace(0.5, 3.5, 11);
	nd = numel(dx);
	for id=1:nd
		tmp = kde_pmf2(x, y, 'dx', dx(id));
		Hv(id) = Hf(tmp); % H for several dx values
	end

	im subplot 3
	plot(dx([1 end]), [h h], '-', ...
		dx, Hv, 'o', ...
		dx, Hv + log(dx.*dx), '+') % H + log(dx*dy) approximately is h
	legend('h', 'H', 'H + log(dx dy)', 'location', 'southwest'), grid
end

if 0
	clf
	while(1)
		im(cent{:}, hist/n)
		axis([-1 1 -1 1]*10)
		prompt
		im(xs, ys, pmf)
		axis([-1 1 -1 1]*10)
		prompt
	end
end
