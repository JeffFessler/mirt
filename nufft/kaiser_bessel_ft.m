 function y = kaiser_bessel_ft(u, J, alpha, kb_m, d)
%function y = kaiser_bessel_ft(u, J, alpha, kb_m, d)
%|
%| Fourier transform of generalized Kaiser-Bessel function,
%|	in dimension d (default 1).
%|	shape parameter "alpha" (default 2.34 J)
%|	order parameter "kb_m" (default 0)
%| See (A3) in lewitt:90:mdi, JOSA-A, Oct. 1990.
%|
%| in
%|	u	[M 1]	frequency arguments
%|			or 'string' to return function string
%|			or 'handle' to return function handle
%|
%| "opt"
%|	J		default 6
%|	alpha		default 2.34 * J
%|	kb_m		default 0
%|	d		default 1
%|
%| out
%|	y	[M 1]	transform values
%|
%| Copyright 2001-3-30, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(u, 'test'), kaiser_bessel_ft_test, return, end

if ~isvar('J'), J = 6; end
if ~isvar('alpha') || isempty('alpha'), alpha = 2.34 * J; end
if ~isvar('kb_m') || isempty('kb_m'), kb_m = 0; end
if ~isvar('d'), d = 1; end

% trick to return functions
if ischar(u)
	kernel_ft = sprintf('kaiser_bessel_ft(t, %d, %g, %g, 1)', ...
		J, alpha, kb_m);

	switch u
	case 'handle'
		y = @(t) kaiser_bessel_ft(t, J, alpha, kb_m, 1);
	case 'string'
		y = kernel_ft;
	case 'inline'
		y = inline(kernel_ft, 't'); % todo: obsolete
	otherwise
		fail('bad argument "%s"', u)
	end
return
end


% Check for validity of FT formula
persistent warned
if (kb_m < -1)
	if isempty(warned) % only print this reminder the first time
		warn('kb_m=%g < -1 in kaiser_bessel_ft()', kb_m)
		warn('validity of FT formula uncertain for kb_m < -1')
		warned = 1;
	end
elseif (kb_m < 0) && ((abs(round(kb_m)-kb_m)) > eps)
	if isempty(warned) % only print this reminder the first time
		warn('Neg NonInt kb_m=%g in kaiser_bessel_ft()', kb_m)
		warn('validity of FT formula uncertain')
		warned = 1;
	end	
end

% trick: simplified since matlab's besselj can handle complex args!
z = sqrt( (2*pi*(J/2)*u).^2 - alpha^2 );
nu = d/2 + kb_m;
y = (2*pi)^(d/2) .* (J/2)^d .* alpha^kb_m ./ besseli(kb_m, alpha) ...
	.* besselj(nu, z) ./ z.^nu;
y = reale(y);

% old inefficient way:
%Lambda = gamma(nu+1) * (z/2).^(-nu) .* besselj(nu, z);
%y2 = (1/2)^kb_m * pi^(d/2) * (J/2)^d * alpha^kb_m ...
%	.* Lambda ./ gamma(d/2+kb_m+1) ./ besseli(kb_m, alpha);


% kaiser_bessel_ft_test
function kaiser_bessel_ft_test
J = 5; alpha = 6.8;
N = 2^10;
x = [-N/2:N/2-1]'/N * (J+3)/2;
dx = x(2) - x(1);
du = 1 / N / dx;
u = [-N/2:N/2-1]' * du;
uu = 1.5*linspace(-1,1,201)';

mlist = [-2 0 2 7];
leg = {};
for ii=1:length(mlist)
	kb_m = mlist(ii);
	yy(:,ii) = kaiser_bessel(x, J, alpha, kb_m);
	Yf(:,ii) = reale(fftshift(fft(fftshift(yy(:,ii))))) * dx;
	Y(:,ii) = kaiser_bessel_ft(u, J, alpha, kb_m, 1);
	Yu(:,ii) = kaiser_bessel_ft(uu, J, alpha, kb_m, 1);
	leg{ii} = sprintf('m=%d', kb_m);
end

if 0
	plot(u, Yf(:,3), 'cx', u, Y(:,3), 'yo', uu, Yu(:,3), 'y-')
	legend('FFT', 'FT coarse', 'FT fine')
	axis tight, axisx(range(uu)), grid
return
end

if im
	plot(	uu, Yu(:,1), 'c-', uu, Yu(:,2), 'y-', ...
		uu, Yu(:,3), 'm-', uu, Yu(:,4), 'g-')
	axis tight, legend(leg)
	hold on, plot(u, Yf(:,2), 'y.'), hold off
	axisx([-1 1]*max(abs(uu)))
	xlabelf '$\nu$'
	ylabelf '$Y(\nu)$'
	titlef('KB FT, \\alpha=%g', alpha)
end
