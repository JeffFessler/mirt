 function [sino, Hk, hn, nn] = fbp2_sino_filter(type, sino, varargin)
%function [sino, Hk, hn, nn] = fbp2_sino_filter(type, sino, [options])
%|
%| Apply ramp-like filters to sinogram(s) for 2D FBP image reconstruction.
%| Both parallel-beam and fan-beam tomographic geometries are supported.
%| This approach of sampling the band-limited ramp avoids the aliasing that
%| would be caused by sampling the ramp directly in the frequency domain.
%|
%| in
%|	type		'arc' (3rd generation CT) or 'flat' (for parallel too)
%|	sino	[nb (L)] sinogram(s)
%| options
%|	dr | ds	(real)	sample spacing (in distance units, e.g., cm) (default 1)
%|	dsd	(real)	source-to-detector distance, for 'arc' case only.
%|	extra		# of extra sinogram radial samples to keep (default: 0)
%|	npad		# of padded samples. (default: 0, means next power of 2)
%|	decon1		deconvolve effect of linear interpolator? (default: 1)
%|	window	[npad]	samples of apodization window function
%|			for [-np/2,...,np/2-1]. (default: '' = plain ramp)
%|			or a string like 'hann' for some predefined windows
%| out
%|	sino	[nb (L)] filtered sinogram rows
%|	Hk	[npad]	apodized ramp filter frequency response
%|	hn	[npad]	samples of band-limited ramp filter
%|	nn	[npad]	[-np/2,...,np/2-1] vector for convenience
%|
%| Copyright 2005-12-19, Jeff Fessler, University of Michigan

if nargin == 1 && streq(type, 'test'), fbp2_sino_filter_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.ds = 1;
arg.dsd = [];
arg.window = '';
arg.extra = 0;
arg.npad = 0;
arg.decon1 = false;
arg = vararg_pair(arg, varargin, 'subs', {'dr', 'ds'; 'Dsd', 'dsd'});

dims = size(sino);
sino = reshape(sino, dims(1), []);
[sino Hk hn nn] = fbp2_sino_filter_do(type, sino, ...
	arg.ds, arg.dsd, arg.window, arg.extra, arg.npad, arg.decon1);
sino = reshape(sino, [size(sino, 1) dims(2:end)]);


%
% fbp2_sino_filter_do()
%
function [sino, Hk, hn, nn] = fbp2_sino_filter_do(type, sino, ...
	ds, dsd, window, extra, npad, decon1);

[nb na] = size(sino);
if ~npad
	npad = 2^ceil(log2(2*nb-1)); % padded size
%	printm('nb=%d npad=%d', nb, npad)
end
sino = [sino; zeros(npad-nb,na)]; % padded sinogram

[hn nn] = fbp_ramp(type, npad, ds, dsd);

Hk = reale(fft(fftshift(hn)));

Hk = Hk .* fbp2_window(npad, window);

Hk = ds * Hk; % differential for discrete-space convolution vs integral

% linear interpolation is like blur with a triangular response,
% so we can compensate for this approximately in frequency domain
if decon1
	Hk = Hk ./ fftshift(nufft_sinc(nn / npad).^2);
end

sino = ifft_sym( fft(sino, [], 1) .* repmat(Hk, [1 na]), [], 1); % apply filter

% trick: possibly keep extra column(s) for zeros!
sino = sino([1:(nb+extra)],:);
sino([(nb+1):(nb+extra)],:) = 0;


%
% test
%
function fbp2_sino_filter_test
nb = 2^5;
sino = zeros(nb,2); sino(nb/2+1,1) = 1;
ds = 1;
dsd = 30;
[sino1 H1 h1 nn] = fbp2_sino_filter('arc', sino, 'ds', ds, 'dsd', dsd);
[sino2 H2 h2 nn] = fbp2_sino_filter('flat', sino, 'ds', ds);
[sino2 H3 h2 nn] = fbp2_sino_filter('flat', sino, 'ds', ds, 'decon1', 1);
%max_percent_diff(h1, h2) % small!
%max_percent_diff(H1, H2) % small!

if im
	clf, subplot(121)
	plot(nn, h1, 'o', nn, h2, '+')
	xlabel 'n', ylabel 'h[n]'
	axis tight
	legend('arc', 'flat')

	subplot(122)
	plot(nn, fftshift(H1), 'o', nn, fftshift(H2), '+', nn, fftshift(H3), '.')
	xlabel 'k', ylabel 'H[k]', axis tight, axisy(0, 0.5)
	legend('arc', 'flat', 'decon', 'location', 'north')
end
