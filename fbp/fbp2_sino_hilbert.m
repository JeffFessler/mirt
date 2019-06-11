  function [sino, hn, nn, Hk] = fbp2_sino_hilbert(sino, varargin)
%|function [sino, hn, nn, Hk] = fbp2_sino_hilbert(sino, [options])
%|
%| Apply band-limited Hilbert-transform filter to 2D sinogram.
%| Frequency response: H(u) = -1i * sign(u) * rect(u/2/umax)
%|
%| in
%|	sino	[nb (L)] sinogram(s)
%|
%| options
%|	dr | ds	(real)	sample spacing (in distance units, e.g., cm) (default 1)
%|	npad		# of padded samples. (default: 0, means next power of 2)
%|	decon1		deconvolve effect of linear interpolator? (default: 0)
%|	window	[npad]	samples of apodization window function
%|			for [-np/2,...,np/2-1]. (default: '' = plain ramp)
%|			or a string like 'hann' for some predefined windows
%| out
%|	sino	[nb (L)] filtered sinogram rows
%|	hn	[npad]	samples of band-limited Hilbert transform filter
%|	nn	[npad]	[-np/2,...,np/2-1] vector for convenience
%|	Hk	[npad]	spectral samples on [0 ... np-1]
%|
%| Copyright 2011-07-16, Jeff Fessler, University of Michigan

if nargin == 1 && streq(sino, 'test'), fbp2_sino_hilbert_test, clear, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.dr = 1;
arg.window = '';
arg.npad = 0;
arg.decon1 = false;
arg = vararg_pair(arg, varargin, 'subs', {'ds', 'dr'});

dims = size(sino);
sino = reshape(sino, dims(1), []);
[sino hn nn Hk] = fbp2_sino_hilbert_do(sino, ...
	arg.dr, arg.window, arg.npad, arg.decon1);
sino = reshape(sino, [size(sino, 1) dims(2:end)]);


% fbp2_sino_hilbert_do()
function [sino, hn, nn, Hk] ...
	= fbp2_sino_hilbert_do(sino, dr, window, npad, decon1);

[nb na] = size(sino);
if ~npad
	npad = 2^ceil(log2(2*nb-1)); % padded size
end
sino = [sino; zeros(npad-nb,na)]; % padded sinogram

[hn nn] = fbp2_filter_hilbert_make(npad, dr);

Hk = 1i * imag_check(fft(fftshift(hn))); % trick: pure imaginary!

Hk = Hk .* fbp2_window(npad, window);

Hk = dr * Hk; % differential for discrete-space convolution vs integral

% linear interpolation is like blur with a triangular response,
% so we can compensate for this approximately in frequency domain
if decon1
	Hk = Hk ./ fftshift(nufft_sinc(nn / npad).^2);
end

sino = ifft_sym( fft(sino, [], 1) .* repmat(Hk, [1 na]), [], 1); % apply filter


% fbp2_filter_hilbert_make()
function [hn, nn] = fbp2_filter_hilbert_make(n, dr)
nn = [-(n/2):(n/2-1)]';
u0 = 1/2/dr;
hn = 2 * u0 * nufft_sinc(nn/2) .* sin(pi/2 * nn);


% imag_check()
% take imaginary part but check that real part is negligible
function out = imag_check(in, varargin)
out = reale(-1i * in, varargin{:});


% fbp2_sino_hilbert_test()
function fbp2_sino_hilbert_test
nb = 2^5;
dr = 0.1;
[sino1 h1 nn H1] = fbp2_sino_hilbert(zeros(nb,2), 'dr', dr);
h2 = 1 ./ (pi * nn * dr); % samples of ideal non-bandlimited

if im
	clf, subplot(121)
	plot(nn, h1, '.-', nn, h2, '-')
	xlabel 'n', ylabel 'h[n]'
	legend('band-limited', 'ideal')

	subplot(122)
	plot(nn, imag_check(fftshift(H1)), '.-', nn, -sign(nn), '-')
	xlabel 'k', ylabel 'imag(H[k])'
	legend('band-limited', 'ideal', 'location', 'north')
end
