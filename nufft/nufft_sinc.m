 function y = nufft_sinc(x)
%function y = nufft_sinc(x)
%|
%| my version of "sinc" function, because matlab's sinc() is in a toolbox
%|
%| Copyright 2001-12-8, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(x, 'test'), nufft_sinc_test, return, end

iz = find(x == 0); % indices of zero arguments
x(iz) = 1;
y = sin(pi*x) ./ (pi*x);
y(iz) = 1;


% test
function nufft_sinc_test

x = linspace(-4, 4, 2^21+1)';

nufft_sinc(0); % warm up
cpu etic
y1 = nufft_sinc(x);
cpu etoc 'nufft_sinc time'

if 2 == exist('sinc')
	sinc(0); % warm up
	cpu etic
	y2 = sinc(x);
	cpu etoc 'matlab sinc time'
	jf_equal(y1, y2)
end

if im, plot(x, y1, '-'), end
