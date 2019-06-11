 function y = ir_fftshift2(x)
%function y = ir_fftshift2(x)
%|
%| 2D fft shift of 2D+T object

if nargin < 1, ir_usage, end

y = fftshift(fftshift(x, 1), 2);
