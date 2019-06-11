 function [Hk, hn, nn] = fbp_filter(type, n, ds, varargin)
%function [Hk, hn, nn] = fbp_filter(type, n, ds, [options])
%
% 'ramp-like' filters for parallel-beam and fan-beam FBP reconstruction.
% This sampled band-limited approach avoids the aliasing that would be
% caused by sampling the ramp directly in the frequency domain.
%
% in
%	type		'arc' (3rd generation CT) or 'flat' (for parallel too)
%	n	(int)	# of samples (must be even)
%	ds	(real)	sample spacing (in distance units, e.g., cm)
% options
%	Dsd	(real)	source-to-detector distance, for 'arc' case
%	window	[n]	or string like 'hann' for some predefined windows
% out
%	Hk	[n]	apodized ramp filter frequency response
%	hn	[n]	samples of band-limited ramp filter
%
% Copyright 2005-12-16, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(type, 'test'), fbp_filter_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.Dsd = [];
arg.window = '';
arg = vararg_pair(arg, varargin);

[hn nn] = fbp_ramp(type, n, ds, arg.Dsd);

Hk = reale(fft(fftshift(hn)));

Hk = fbp_apodize(Hk, n, arg.window);


%
% fbp_apodize()
%
function H = fbp_apodize(H, n, window)

if ischar(window)
	if isempty(window) || streq(window, 'ramp')
		window = ones(n,1);
	elseif streq(window, 'hann')
		window = hann(n, 'periodic');
	elseif streq(window, 'hann50')
		window = zeros(n,1);
		window(n/4+1:3*n/4) = hann(n/2, 'periodic');
	elseif streq(window, 'hann75')
		window = zeros(n,1);
		window(n/8+1:7*n/8) = hann(3*n/4, 'periodic');
	else
		error 'unknown window'
	end
elseif length(window) ~= n
	error 'bad window length'
end

H = H .* fftshift(window);

%
% test
%
function fbp_filter_test
nb = 2^5;
ds = 1;
Dsd = 10;
[H1 h1 nn] = fbp_filter('arc', nb, ds, 'Dsd', Dsd);
[H2 h2 nn] = fbp_filter('flat', nb, ds);
max_percent_diff(h1, h2)
max_percent_diff(H1, H2)

clf, subplot(121)
plot(nn, h1, 'o', nn, h2, '+')
xlabel 'n', ylabel 'h[n]'
axisx(-10, 10)
legend('arc', 'flat')

subplot(122)
plot(nn, fftshift(H1), 'o', nn, fftshift(H2), '+')
xlabel 'k', ylabel 'H[k]'
%axisx(-10, 10)
legend('arc', 'flat')
