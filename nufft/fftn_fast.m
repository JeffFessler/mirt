 function Xs = fftn_fast(xs, ns)
%function Xs = fftn_fast(xs, ns)
%|
%| For some reason, in matlab versions before about 7.4 (R2007a),
%| matlab's fftn routine was suboptimal for the case of 2D FFTs,
%| at least on some machines.
%| The improvement herein was found by Hugo Shi.
%| After 7.4, fftn worked fine so this routine is no longer needed.
%|
%| Note: At ISMRM 2009, Phil Beatty mentioned that sequential 1D FFT approach
%| can be faster particularly when "trimming" of excess FOV is used.
%|
%| Copyright 2004-6-28, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), return, end
if streq(xs, 'test'), fftn_fast_test, return, end

% around version 7.4, fftn was fastest, so hardwire that!
if nargin < 2
	Xs = fftn(xs);
else
	if length(ns) == 1
		Xs = fft(xs, ns);
	else
		Xs = fftn(xs, ns);
	end
end
return


if nargin < 2, ns = size(xs); end

if ndims(xs) == 2 % 2D or 1D case
	if min(size(xs)) == 1
		if ns(2) ~= 1, error 'bug', end
		Xs = fft(xs, ns(1));
	else
		Xs = fftn_fast_fftfft(xs, ns);
	end
else
	Xs = fftn(xs, ns);
end


function Xs = fftn_fast_fftfft(xs, ns)
Xs = fft(fft(xs, ns(1)).', ns(2)).';


% test configuration of fftn_fast for this machine
function fftn_fast_test
fftn_fast_test2 % 2D


% test configuration of fftn_fast for this machine for 2D
function fftn_fast_test2

n = 2^8;
x = rand(n,n);
ns = [2*n n];

printf('starting test; be patient.')

% first loop is to get everything in cache or whatever.
% doing it twice is the only way to get an accurate comparison!
for nloop = [2 40];
	tic, for ii=1:nloop, X{1} = fftn_fast(x, ns); end
	tt(1) = toc; ty{1} = 'fftn_fast';

	tic, for ii=1:nloop, X{2} = fftn(x, ns); end
	tt(2) = toc; ty{2} = 'fftn';

	tic, for ii=1:nloop, X{3} = fft(fft(x, ns(1)).', ns(2)).'; end
	tt(3) = toc; ty{3} = 'fftfft_inline';

	tic, for ii=1:nloop, X{4} = fft(fft(x, ns(1), 1), ns(2), 2); end
	tt(4) = toc; ty{4} = 'fftfft_brack';

	tic, for ii=1:nloop, X{5} = fft2(x, ns(1), ns(2)); end
	tt(5) = toc; ty{5} = 'fft2';

	tic, for ii=1:nloop, X{6} = fftn_fast_fftfft(x, ns); end
	tt(6) = toc; ty{6} = 'fftfft_func';
end

for ii = 1:length(tt)
	printf('time %14s = %g', ty{ii}, tt(ii))
	if max_percent_diff(X{1}, X{ii}) > 1e-11, error 'bug', end
end

if tt(1) > 1.05 * min(tt(2:end))
	warn 'fftn_fast may be configured supoptimally for your machine!'
	printf('fftn / fftn_fast = %g%% ', tt(2) / tt(1) * 100.)
else
	printf('fftn_fast is configured appropriately for your machine')
	printf('fftn / fftn_fast = %g%% ', tt(2) / tt(1) * 100.)
end
