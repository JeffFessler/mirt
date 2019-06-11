 function w = ir_hann_periodic(M)
%function w = ir_hann_periodic(M)
%|
%| Hanning window with 0 only in the first element.
%| Useful with FFT
%| Equivalent to matlab hann(M, 'periodic')
%| Needed because octave 3.6.4 does not offer periodic option.
%|
%| 2013-08-10, Jeff Fessler, University of Michigan

w = 0.5 * (1 - cos(2*pi*(0:M-1)'/M));
