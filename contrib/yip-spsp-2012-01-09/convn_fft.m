function y = convn_fft(x,f)
%function y = convn_fft(x,f)
%this functional performs y = convn(x,f), but is implemented in the
%frequency domain. x is the input signal or image, and f is the kernel. The
%inputs can be either 2D or 3D. 
%cy 11/3/2008
%

disp('Filtering via N-dimensional FFT...');
X = fftshift(fftn(fftshift(x)));
F = abs(fftshift(fftn(f,size(X))));
y = fftshift(ifftn(fftshift(F.*X)));

