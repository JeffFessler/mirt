  function [recon_ph_demod, recon_im, lp_im, full_kspace] = ...
	homodyne_recon(partial_kspace, N, M, overlap, varargin)
%|function [recon_im, recon_ph_demod, lp_im, full_kspace] =
%|	homodyne_recon(partial_kspace, N, M, overlap, varargin)
%|
%| Reconstruct image from partial kspace data via homodyne method outlined
%| in Doug Noll's 1991 paper "Homodyne Detection in Magnetic Resonance Imaging"
%|
%| in
%|	partial_kspace	[N/2+overlap+1 M] for direction=1, or
%|			[N M/2+overlap+1] for direction=2
%|			2D DFT of original image, with almost half zeroed out
%|			(additional strip of samples of width overlap beyond DC)
%|	N		dimension 1 of full kspace
%|	M		dimension 2 of full kspace
%|	overlap		width of additional strip of samples beyond DC,
%|			used in low-passed image for phase correction
%| option:
%|	'direction'	1 (default) or 2
%|	'merge'		'ramp' (default) or 'step' (merging filter)
%|
%| out
%|	recon_ph_demod	[N M]	reconstructed image, demodulated by reference
%|
%| out (for debugging only)
%|	recon_im	[N M]	reconstructed from partial kspace data,
%|			(does not use low pass reference image to demodulate)
%|	lp_im		[N M]	low-pass reference image
%|				(used to demodulate reconstructed image)
%|	full_kspace	[N M]	reconstructed full kspace
%|
%| 2012-06-16 Mai Le, University of Michigan
%| 2012-10-15 tweaks by JF

if nargin == 1 && streq(partial_kspace, 'test')
	run_mfile_local('homodyne_recon_demo')
return
end
if nargin < 4, help(mfilename), error(mfilename), end

arg.direction = 1;
arg.merge = 'ramp';
arg = vararg_pair(arg, varargin);

if arg.direction == 2 % trick: switch dimensions
	[recon_ph_demod, recon_im, lp_im, full_kspace] = ...
		homodyne_recon(partial_kspace.', M, N, overlap, ...
			'direction', 1, 'merge', arg.merge);
	recon_ph_demod = recon_ph_demod.';
	recon_im = recon_im.';
	lp_im = lp_im.';
	full_kspace = full_kspace.';
return

elseif arg.direction ~= 1
	fail('bad direction %d', arg.direction)
end

if (N < size(partial_kspace,1)) || (M ~= size(partial_kspace,2))
	fail('[N M]=[%d %d] but size(partial_kspace)=[%d %d]', ...
		N, M, size(partial_kspace,1), size(partial_kspace,2))
end

tmp = zeros(N-size(partial_kspace,1),M);
partial_kspace = [partial_kspace; tmp];

% construct merging filter, can use either ramp or step
switch arg.merge
case 'ramp'
	mask = homodyne_recon_ramp(N, overlap);
	mask = repmat(mask', 1, M);
case 'step'
	mask = homodyne_recon_step(N, overlap);
	mask = repmat(mask', 1, M);
otherwise
	fail('bad merge %s', arg.merge)
end
%im(1, mask), im(2, partial_kspace), prompt

full_kspace = mask .* partial_kspace;
recon_im = ifft2(ifftshift(full_kspace));
lpf = get_lpf(N, M, overlap);
lp_im = ifft2(ifftshift(lpf .* partial_kspace));
ph_hat = angle(lp_im);
interm = recon_im .* exp(-1i*ph_hat);
recon_ph_demod = real(interm);

end % homodyne_recon()


% homodyne_recon_ramp()
function ramp = homodyne_recon_ramp(N, overlap)
if (mod(N,2) == 0) % for N even
	ramp = [1 2*ones(1,N/2-overlap-1) linspace(2,0,min(overlap*2+1,N-1)) ...
		zeros(1,N/2-overlap-1)];
else
	ramp = [2*ones(1,(N-1)/2-overlap) linspace(2,0,min(overlap*2+1,N)) ...
		zeros(1,(N-1)/2-overlap)];
end
end % homodyne_recon_ramp()


% homodyne_recon_step()
function step = homodyne_recon_step(N, overlap)
if (mod(N,2) == 0) % for N even
	step = [1 2*ones(1,N/2-overlap-1) ones(1,min(overlap*2+1,N-1)) ...
		zeros(1,N/2-overlap-1)];
else
	step = [2*ones(1,(N-1)/2-overlap) ones(1,min(overlap*2+1,N-1)) ...
		zeros(1,(N-1)/2-overlap)];
end
end % homodyne_recon_step()


% get_lpf()
% Generate softer low pass filter to reduce ringing
function lpf = get_lpf(N, M, overlap)
if (mod(N,2) == 0)
	d1 = max(N/2+1-overlap,1):min(N/2+1+overlap,N);
else
	d1 = max((N+1)/2-overlap,1):min((N+1)/2+overlap,N);
end
d2 = 1:M;
sigma = length(d1)/5;
lpf = zeros(N,M);
%lpf(d1,d2) = fspecial('gaussian', [length(d1) length(d2)], sigma); % in IP toolbox
lpf(d1,d2) = my_gauss([length(d1) length(d2)], sigma);
lpf = lpf ./ max(lpf(:)); % normalize

% todo: provide option to apodize only in the truncated direction!

end % get_lpf()


% my_gauss()
% generates 2D Gaussian filter of given dims and st dev
function gfilter = my_gauss(dims, sigma)
[xx,yy] = ndgrid(1:dims(1),1:dims(2));
% center values
xx = xx - dims(1)/2;
yy = yy - dims(2)/2;
gfilter = exp(-(xx.^2+yy.^2)/(2*sigma^2))/(2*pi*sigma^2);
end % my_gauss
