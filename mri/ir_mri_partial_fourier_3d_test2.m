% ir_mri_partial_fourier_3d_test2.m
% test 2D partial Fourier MRI method ir_mri_partial_fourier_3d.m
% to examine behavior of xu:01:pfi
% 2014-06-? Mai Le
% 2014-06-22 JF mods

%% --- Parameters that control recon quality

percentage = 0.65; % partial k-space fraction, must be greater than 0.5

window_step3 = 12; % controls transition band width for apodization
window_step8 = 6;


%% --- Chose test case

nz = 1;
% uncomment line below to test 3D
% even number just big enough to pass condition in PF_3D that requires apodization area of 16:
%nz = 2*ceil((2*max(window_step3,window_step8)+1) / (2*(percentage - 0.5)) / 2);

% Shepp-Logan
%xtrue = phantom('Modified Shepp-Logan');
if nz > 1
	xtrue = ellipse_im(64); % need even dims
	xtrue = xtrue(5:64,:); % stress test non-square
	xtrue = repmat(xtrue, [1 1 nz]);
im(xtrue) % fails in octave - todo try in matlab
return
	xtrue(:,:,1) = zeros(size(xtrue(:,:,1))); xtrue(10, 50) = 1; % impulse
else
	xtrue = ellipse_im(128);
	xtrue = xtrue(5:124,:); % stress test non-square
%	xtrue = zeros(size(xtrue)); xtrue(10, 80) = 1; % impulse
end
[nx ny nz] = size(xtrue)

% make it complex
xtrue = xtrue .* exp(-1i*0.05); % constant phase - simple case
xtrue = xtrue .* exp(-1i*2*pi*(([1:nx]'*[1:ny])/(nx*ny)).^2); % nonlinear phase

% Tough Case
%xtrue = 1 + rand(nx, ny, nz) .* exp(1i*0.05); % not quite constant phase :(

kspace = fftshift(fftn(xtrue));

if 0 % Artificial k-space for visualizing effect of algo
	kspace = rand(size(xtrue));
	xtrue = ifftn(ifftshift(k));
end


%% --- Select partial k-space and do PF recon

pf_location = [1 0 0]; % stress test
nx_pf = round(percentage*nx);
ny_pf = round(percentage*ny);
nz_pf = round(percentage*nz);
kx = (1:nx_pf) + pf_location(1) * (nx - nx_pf);
ky = (1:ny_pf) + pf_location(2) * (ny - ny_pf);
kz = (1:nz_pf) + pf_location(3) * (nz - nz_pf);

sampling = false(size(xtrue));
sampling(kx, ky, kz) = true;
% partial_kspace = sampling .* k;
partial_kspace = kspace(kx, ky, kz);

% log scale over 6 digits (single precision)
im_log = @(i, x, t) ...
	im(i, log10(abs(x) / max(abs(kspace(:))) + 1e-20), [-6 0], t);

ir_fontsize im_axes 5
ir_fontsize title 10
ir_fontsize label 10
im plc 3 5
im(1, abs(xtrue), 'true')
im(2, sampling)
tmp = @(k, n) unique([1 n minmax(k)']);
xtick(tmp(kx, nx))
ytick(tmp(ky, ny))
im_log(6, kspace, 'full k-space')
im_log(7, embed(partial_kspace(:), sampling), 'partial')

if nz > 1
	full_dims = [nx ny nz];
else
	full_dims = [nx ny];
end

if ~isvar('img_est')
	[img_est, full_kspace] = ir_mri_partial_fourier_3d(...
		partial_kspace, full_dims, 'pf_location', pf_location, ...
		'chat', 1, 'show', 1, 'niter', 10, ...
		'window_step3', window_step3, 'window_step8', window_step8);
%		'init', xtrue, ... % of course it works
end

im(5, abs(img_est)), title '|img_est|'
im(10, angle(img_est)), title '\angle < img_est'

diff_img = img_est - xtrue;
%im(9, abs(img), '|Recon|'), cbar
im(11, abs(diff_img), 'Error')
xlabelf('max = %.3g', max(abs(diff_img(:))))
im_log(15, full_kspace, 'kspace est')
im_log(12, full_kspace - kspace, 'kspace err')

if 0
	im plc 1 1
	im_log(1, full_kspace - kspace, 'kspace err')
end

%im_log(11, full_kspace, 'estimated')

% conclusions:
% - perfect recovery not possible
% - missing corners of k-space filled in by convolution of fft of phase difference between current iterate and low frequency estimate
