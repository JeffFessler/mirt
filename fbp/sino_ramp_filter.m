 function sino = sino_ramp_filter(sino)
%function sino = sino_ramp_filter(sino)
%|
%| apply ramp filter to a sinogram using simple FFT method
%|
%| use fbp2_sino_filter.m instead of this!
%|
%| in
%|	sino	[nb na] sinogram
%| out
%|	sino	[nb na]	ramp filtered sinogram
%|
%| Copyright 2001-10-4, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

[nb na] = size(sino);
npad = 2*nb;
if 2^round(log(nb)/log(2)) ~= nb
	warning 'fix: pad for power-of-2 for efficiency'
end

sino = [sino; zeros(size(sino))]; % zero pad
ramp = abs([-nb:nb-1]') / (2*nb); % pure ramp filter
ramp(ramp == 0) = 2/pi^2 / (2*nb); % DC fix from crawford:91:cfa
sino = ifft(fft(sino) .* fftshift(repmat(ramp, [1 na])));
sino = real(sino);
sino = sino(1:nb,:);
