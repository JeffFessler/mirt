 function y = ifft_sym(varargin)
%function y = ifft_sym(varargin)
%|
%| matlab 7.0 introduced a 'symmetric' option to ifft to handle
%| spectra that are (circularly) hermitian symmetric (real signal).
%| this glue routine is to provide backward compatibility for matlab 6.5.
%| Caution: v7 ifft with 'symmetric' just uses the first half of the spectrum
%| along whichever dimension is requested.  Here, for pre v7, I just take
%| the real part.  The difference is neglible in the cases where this
%| routine is expected to be used, where the spectrum should be exactly
%| symmetric but has slight asymmetry due to numerical precision.
%| If the spectrum is severely asymmetric, then "real(ifft())" and
%| ifft(..., 'symmetric') will differ substantially.  (But one should
%| not call this routine in such cases.)

if ~nargin, ir_usage, end
if nargin == 1 && streq(varargin{1}, 'test'), ifft_sym_test, return, end

if ir_is_octave || is_pre_v7
	y = ifft(varargin{:});
	if isa(varargin{1}, 'double')
		tol = 1e-11;
	else
		tol = 1e-6;
	end

	y = reale(y, tol, 'prompt');
else
	y = ifft(varargin{:}, 'symmetric');
end

function y = ifft_sym_test
del = 10^5*eps;
format compact
x1 = [4 2+0i*del 8 2-1i*del]
y1 = ifft(x1)
y2 = ifft_sym(x1)
x2 = fft(y2)
y1 - y2
