 function t = min_cos_quad(m, p, b, c, niter)
%function t = min_cos_quad(m, p, b, c, niter)
%|
%| find the minimizer of the sinusoid + quadratic form:
%| f(t) = m * (1-cos(t-p)) + b*t + c/2 t^2
%|
%| This is useful for certain phase-related optimization transfer problems.
%|
%| Copyright 2003-11-9, Jeff Fessler, University of Michigan

if nargin < 4
	if nargin == 1 && streq(m, 'test')
		min_cos_quad_test, return
	elseif nargin == 1 && streq(m, 'fig')
		min_cos_quad_fig, return
	else
		ir_usage
	end
end

if any(c < 0), fail 'need c > 0', end
if any(m < 0), fail 'need m >= 0', end

t = zeros(size(m));

% analytical solution where m=0
i0 = m == 0;
t(i0) = -b(i0) ./ c(i0);

p0 = p(~i0);
a = -p0 - b(~i0) ./ c(~i0);
k = c(~i0) ./ m(~i0);

u = zero_sin_line(-p0, k, a, niter);
t(~i0) = u + p0;


% zero_sin_line()
% find zeros of sin(u) + k * (u - a)
function u = zero_sin_line(u, k, a, niter)

amid = round(a/(2*pi)) * 2*pi; % midpoint of nearest [-pi,pi] + n*2*pi
% for any initial points not within proper interval,
% set initial guess to that interval's midpoint
ibad = (u < amid - pi) | (u > amid + pi);
u(ibad) = amid(ibad);

for ii=1:niter
%	s = mod(t+pi,2*pi) - pi; % [-pi,pi]
%	curv = nufft_sinc(s/pi);
%	denom = m .* curv + c;
%	t = t - (m .* sin(t-p) + b + c.*t) ./ denom;

%	s = mod(t+pi,2*pi) - pi; % [-pi,pi]
%	curv = nufft_sinc(s/pi);
%	denom = curv + k;
	denom = 1 + k;
	u = u - (sin(u) + k.*(u-a)) ./ denom;
end


% min_cos_quad_test
% self test
function min_cos_quad_test

mlist = 1;
plist = 3.8;
blist = 0.1;
clist = 0.2;
[mm pp bb cc] = ndgrid(mlist, plist, blist, clist);
im clf, pl=100 + 10*numel(mm);
t = linspace(-3*pi,3*pi,501);
f = @(t, m, p, b, c) m * (1-cos(t-p)) + b*t + c/2 * t.^2;
for ip=1:numel(mm)
	m = mm(ip);
	p = pp(ip);
	b = bb(ip);
	c = cc(ip);

	if im
		subplot(pl+ip)
		plot(t, f(t,m,p,b,c), 'y-', 0, f(0,m,p,b,c), 'go')
		xlabel 't', ylabel 'surrogate(t)'
		title 'm [1-cos(t-p)] + b t + c/2 t^2'
		axis([minmax(t)' 0 9])
		xaxis_pi '-3p -p 0 p 3p'
	end
	for ii=1:4
		x = min_cos_quad(m, p, b, c, ii);
		if im
		hold on, plot(x, f(x,m,p,b,c), 'co'), hold off
		end
	end
	if im
		hold on, plot(x, f(x,m,p,b,c), 'r.'), hold off
	end
	% ir_savefig c 'fig_cos_quad'
end


%
% figure to show minima
%
function min_cos_quad_fig

alist = [pi/4 3*pi/4 -pi/2]+0*pi;
klist = [0.1 0.2 0.2];
t = linspace(-2*pi,2*pi,301);
d1 = @(x, a, b, c) a * (1-cos(x)) + c/2 * (x-b).^2;
clf, pl = @(i) subplot(330 + i);
for ii=1:3
	a = alist(ii);
	k = klist(ii);

	pl(ii), plot(t, -cos(t)+k/2*(t-a).^2, '-')
	xaxis_pi '-2p -p 0 p 2p'

	pl(ii+3), plot(t, -sin(t), '--', t, k*(t-a), '--')
	xaxis_pi '-2p -p 0 p 2p'
	axisy(-2,2), ytick([-2 0 2])
	if ii==3, legend('-sin(t)', 'k(t-a)',2), end

	pl(ii+6), plot(t, cos(t)+k, '-', t, 0*t, ':')
	xaxis_pi '-2p -p 0 p 2p'
	axisy(-2,2), ytick([-2 0 2])
end
