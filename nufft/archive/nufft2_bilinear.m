 function X = nufft2_bilinear(omega, x, K1, K2, n_shift)
%function X = nufft2_bilinear(omega, x, K1, K2, n_shift)
%	in:
%		omega	[M 2]		frequencies in radians
%		x	[N1 N2]		image dimensions
%		J			# of neighbors used
%		K1,K2			FFT sizes (should be > N1,N2)
%		n_shift [2]		n = 0-n_shift to N-1-n_shift
%	out:
%		X	[M 1]	DTFT2 at omega, approximated by bilinear interp
%
%	Like fft(), this expects the signals to be x(0,0), ...
%	Use n_shift = [N1/2 N2/2] for x(-N1/2,-N2/2), ...
%
% Copyright 2001-9-20, Jeff Fessler, University of Michigan

warn('This function is obsolete.  Use nufft_init with the "linear" option')

% if no arguments, then run a simple test
if nargin < 1
	N1 = 4;
	N2 = 8;
	n_shift = [2.7 3.3];
	x = [[1:N1]'*ones(1,3), ones(N1,N2-3)]; % test signal
%	x = zeros(N1,N2); x(1,1) = 1;
	if 0	% test with uniform frequency locations
		o1 = 2 * pi * [0:(N1-1)]' / N1;
		o2 = 2 * pi * [0:(N2-1)]' / N2;
		[o1, o2] = ndgrid(o1, o2);
		omega = [o1(:) o2(:)];
	else	% nonuniform frequencies
		o1 = [0 7.2 2.6 3.3]';
		o2 = [0 4.2 -1 5.5]';
		omega = [o1(:) o2(:)];
	end
	Xd = dtft2(x, omega, n_shift);
	Xb = nufft2_bilinear(omega, x, 2*N1, 2*N2, n_shift);
%	Xb = reshape(Xb, N1, N2)
%	disp([Xd Xb Xb-Xd])
	help(mfilename)
	disp(sprintf('max %% difference = %g', max_percent_diff(Xd,Xb)))
	return
end


if ~isvar('n_shift') || isempty(n_shift), n_shift = [0 0]; end
if ~isvar('useloop') || isempty(useloop), useloop = 0; end

M = size(omega,1);
[N1 N2] = size(x);

koff1 = bilinear_offset(omega(:,1), K1);	% [M 1] nearest
koff2 = bilinear_offset(omega(:,2), K2);	% [M 1] nearest

% indices into oversampled FFT components
k1 = mod(outer_sum(koff1, [0 1]), K1) + 1;	% [M 2] {1,...,K1}
k2 = mod(outer_sum(koff2, [0 1]), K2) + 1;	% [M 2] {1,...,K2}

% 1D interpolation coefficient vectors
[u1l u1r] = interp_coef(omega(:,1), koff1, K1, N1);	% [M 1]
[u2l u2r] = interp_coef(omega(:,2), koff2, K2, N2);	% [M 1]


% precorrect for bilinear interpolation filtering
if 0
	warning 'applying precorrection'
	[i1, i2] = ndgrid([0:(N1-1)]-n_shift(1), [0:(N2-1)]-n_shift(2));
	tmp = (nufft_sinc(i1 / K1) .* nufft_sinc(i2 / K2)).^2;
	printm('sinc correction min=%g', min(tmp(:)))
	x = x ./ tmp;
end

% FFT and bilinear interpolation
% can't use interp2 here because we want periodic end conditions!
% oh, maybe i could with a little padding of ends... 
Xk = fft2(x, K1, K2);
kk11 = k1(:,1) + (k2(:,1) - 1) * K1;
kk21 = k1(:,2) + (k2(:,1) - 1) * K1;
kk12 = k1(:,1) + (k2(:,2) - 1) * K1;
kk22 = k1(:,2) + (k2(:,2) - 1) * K1;
t1 = u1l .* Xk(kk11) + u1r .* Xk(kk21);
t2 = u1l .* Xk(kk12) + u1r .* Xk(kk22);
X = u2l .* t1 + u2r .* t2;

% apply phase shift
phase = exp(i * (omega * n_shift(:)));	% [M 1]
X = X .* phase;


% index from 0 to K-1 (will modulo later)
function koff = bilinear_offset(om, K)
gam = 2*pi/K;
koff = floor(om / gam);

% make 1D interpolation coefficient vectors
function [ul, ur] = interp_coef(om, koff, K, N)
gam = 2*pi/K;
ul = 1 - (om/gam - koff);
ur = 1 - ul;
ul = ul .* kphase(om - koff * gam, N);
ur = ur .* kphase(om - (koff+1) * gam, N);

% natural phase function. trick: force it to be 2pi periodic
function ph = kphase(om, N)
ph = exp(-i * mod0(om,2*pi) * (N-1)/2);
