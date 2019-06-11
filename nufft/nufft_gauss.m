 function [kernel, kernel_ft] = nufft_gauss(type, J, sig)
%function [kernel, kernel_ft] = nufft_gauss(type, J, sig)
%|
%| Gaussian bell kernel functions truncated to support [-J/2,J/2]
%| for NUFFT interplation, with width parameter "sig"
%| in
%|	type		'string' or 'inline'
%|	J		interval width
%|	sig		width parameter
%| out
%|	kernel(k,J)
%|	kernel_ft(t)
%|
%| Copyright 2002-7-15, Jeff Fessler, University of Michigan

% if no arguments, give an example, comparing ft to zn
if nargin < 1
	help(mfilename)
	N = 256; K = 2*N; J = 4;
	[kernel, kernel_ft] = nufft_gauss('inline', J);

	n = [0:N-1]' - (N-1)/2;		% trick due to complex phase term
	sn_ft = 1 ./ kernel_ft(n/K);
	sn_zn = reale(1 ./ nufft_interp_zn(0, N, J, K, kernel));

	k = linspace(-J/2-1,J/2+1,101);
	clf, subplot(121), plot(k, kernel(k,J)), axis tight
	xlabel 'k', ylabel '\psi(k)', title 'Gaussian bell'
	subplot(122), plot(n, sn_ft, 'c-o', n, sn_zn, 'y-')
	axis tight
	legend('1/FT', '1/zn')
	xlabel 't', ylabel '\Psi(t)', title 'Reciprocal of Fourier transform'

	error(mfilename)
end

if ~isvar('type') || isempty(type), type = 'string'; end
if ~isvar('J'), J = 6; end
if ~isvar('sig') || isempty('sig'), sig = 0.78 * sqrt(J); end

if ~ischar(type), error 'type is string', end

kernel = sprintf('exp(-(k/%g).^2/2) .* (abs(k) < J/2)', sig);
kernel_ft = sprintf('%g*sqrt(2*pi)*exp(-pi*(t*%g*sqrt(2*pi)).^2)', sig, sig);


if streq(type, 'string')
	return
elseif streq(type, 'inline')
	kernel = inline(kernel, 'k', 'J');
	kernel_ft = inline(kernel_ft, 't');
else
	error 'type must be "inline" or "string"'
end
