 function [proj, J, alpha, kb_m, d] = kaiser_bessel_xray(r, J, alpha, kb_m, d)
%function [proj, J, alpha, kb_m, d] = kaiser_bessel_xray(r, J, alpha, kb_m, d)
%
% X-ray transform of generalized Kaiser-Bessel function,
% See (A7) in lewitt:90:mdi, JOSA-A, Oct. 1990.
%
% in
%	r	[?]	radial locations in projection space (unitless)
%
% options
%	J		diameter of blob (a = J/2), default 4
%	alpha		shape parameter, default 10.83
%	kb_m		order parameter, default 2
%	d		dimension, default 2
%
% out
%	proj	[?]	x-ray transform values
%
% Copyright 2005-7-29, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(r, 'test'), kaiser_bessel_xray_test, return, end

if ~isvar('J'), J = 4; end
if ~isvar('alpha') || isempty('alpha'), alpha = 10.83; end
if ~isvar('kb_m') || isempty('kb_m'), kb_m = 2; end
if ~isvar('d'), d = 2; end
tol = 1e-4;

%
% trick to yield inline functions
%
if ischar(r)
	kernel = sprintf('kaiser_bessel_xray(r, %d, %g, %g, %g)', ...
		J, alpha, kb_m, d);

	if streq(r, 'string')
		proj = kernel;
	elseif streq(r, 'inline')
		proj = inline(kernel, 'r');
	else
		error 'bad argument'
	end
return
end


%
% Check for validity of FT formula
%
persistent warned
if (kb_m < 0 || (abs(round(kb_m)-kb_m) > eps))
	if isempty(warned)
		printf([mfilename: 'kb_m=%g in kaiser_bessel_xray()'], kb_m)
		printf('validity of formula uncertain')
		warned = 1;
	end
end

a = J/2;
factor = a / besseli(kb_m, alpha) * sqrt(2*pi/alpha);
root = sqrt(1 - (r/a).^2);
nu = kb_m + 1/2;
proj = factor * root.^nu .* besseli(nu, alpha * root);
proj = reale(proj, tol);


%
% test
%
function kaiser_bessel_xray_test
x = linspace(-2.5,2.5,501)';
y = x;
[xx yy] = ndgrid(x, y);
r = sqrt(xx.^2 + yy.^2);
[proj J kb_m alpha] = kaiser_bessel_xray(x);
dx = x(2) - x(1);

f0 = kaiser_bessel(r, J, kb_m, alpha);
psum = sum(f0,2) * dx;

if im
	clf
	im(221, x, y, f0, 'blob'), grid
	subplot(222)
	plot(x, f0(:,y==0)), title 'profile', axis tight
	xtick([-2:1.0:2]), grid
	subplot(212)
	plot(x, proj, '-', x, psum, '--'), axis tight, title 'projection'
	legend('Analytical', 'Numerical')
end
max_percent_diff(proj, psum, mfilename)
