 function y = ir_sinc_nopi(x)
%function y = ir_sinc_nopi(x)
%|
%| no pi version of "sinc" function, because matlab's sinc() is in a toolbox
%|
%| Copyright 2001-12-8, Jeff Fessler, University of Michigan
%| Modified by M Allison to not have a pi.

if nargin < 1, help(mfilename), error(mfilename), end
if streq(x, 'test'), ir_sinc_nopi_test, return, end

iz = find(x == 0); % indices of zero arguments
x(iz) = 1; % arbitrary value
y = sin(x) ./ x;
y(iz) = 1;


% test
function ir_sinc_nopi_test

x = linspace(-4, 4, 2^21+1)' * pi;

ir_sinc_nopi(0); % warm up
cpu etic
y1 = ir_sinc_nopi(x);
cpu etoc 'ir_sinc_nopi time'

if 2 == exist('sinc')
	sinc(0); % warm up
	tmp = x / pi;
	cpu etic
	y2 = sinc(tmp);
	cpu etoc 'matlab sinc time'
	equivs(y1, y2)
end

if im, plot(x/pi, y1, '-'), end
