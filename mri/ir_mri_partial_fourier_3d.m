 function [img_est, full_kspace] = ...
	ir_mri_partial_fourier_3d(partial_kspace, full_dims, varargin)
%function [img, full_kspace] = ...
%|	ir_mri_partial_fourier_3d(partial_kspace, full_dims, varargin)
%|
%| 2D and 3D versions of Xu and Haacke's 2001 "Partial Fourier imaging in
%| multi-dimensions: A means to save a full factor of two in time", JMRI
%| (xu:01:pfi) doi 10.1002/jmri.1228
%|
%| in:
%|	partial_kspace: [pf_nx pf_ny pf_nz]
%|		assumed taken from origin-centered k-space
%|	full_dims [nx ny nz] : dims of full_kspace
%|
%| varargin:
%|	pf_location:
%|		3x1 logical, placement of partial_kspace in full_kspace
%|		0: 1:N
%|		1: (end-N+1):end
%|		default: [1 0 1], for dynamic phantom data 2014/05/08
%|	niter: # iterations of POCS method; default = 5
%|	window_step3
%|		width of transition band for apodization in low freq estimation
%|		default = 8 (set to 0 to disable tapering)
%|	window_step8
%|		width of transition band for apodization in iterative step
%|		default = 3 (set to 0 to disable tapering)
%|	fill_conj 0|1
%|		initialize other "half" of k-space with conjugate (default 0)
%|	init	provide optional initial image (default [])
%|	show 0|1
%|		visualize intermediate steps? default 0
%|
%| out:
%|	img: [nx ny nz]
%|	full_kspace: [nx ny nz]
%|
%| 2014-06-17 Mai Le
%| 2014-06-22 tweaks by JF to work with octave
%|
%| notes:
%|	- "filter" is more of an apodization, reduces discontinuities in k-space
%|	- handles only even sizes for now

if nargin == 1 && streq(partial_kspace, 'test')
	run_mfile_local ir_mri_partial_fourier_3d_test2 % 2D test
	return
end
if nargin < 2, help(mfilename), error(mfilename), end

% options
arg.pf_location = [1 0 1];
arg.niter = 5;
arg.window_step3 = 8;
arg.window_step8 = 3;
arg.fill_conj = false; % initialize other half of k-space with complex conj ?
arg.chat = false;
arg.show = false;
arg.init = [];
arg = vararg_pair(arg, varargin);

% log scale over 6 digits (single precision)
im_log = @(i, x) ...
	im(i, log10(abs(x) / max(abs(x(:))) + 1e-20), [-6 0]);

% detect 2D or 3D
assert(all(mod(full_dims,2) == 0), 'err: odd dims for full kspace')
assert(ndims(partial_kspace) == numel(full_dims), ...
	'full_dims does not match dims of partial_kspace')
switch numel(full_dims)
case 2
	is_3d = 0;
case 3
	is_3d = 1;
otherwise
	fail('invalid number of dimensions: must be 2D or 3D')
end

% determine sizes of full, partial, and low-frequency kspace
nx = full_dims(1);
ny = full_dims(2);
[pf_nx, pf_ny, pf_nz] = size(partial_kspace);
if nx == pf_nx
	lf_kx = 1:pf_nx
else
	lf_kx = ((nx - pf_nx) + 2):pf_nx; % trick: odd to ensure Herm. sym.
end
if ny == pf_ny
	lf_ky = 1:pf_ny;
else
	lf_ky = ((ny - pf_ny) + 2):pf_ny; % trick: odd to ensure Herm. sym.
end
%minmax(lf_kx - (nx/2+1)) % check
if is_3d
	nz = full_dims(3);
	if nz == pf_nz
		lf_kz = 1:pf_nz;
	else
		lf_kz = ((nz - pf_nz) + 2):pf_nz; % trick: odd for Herm. sym.
	end
else
	nz = 1;
	lf_kz = 1;
end

% indeces for partial Fourier
pf_kx = (1:pf_nx) + arg.pf_location(1) * (nx - pf_nx);
pf_ky = (1:pf_ny) + arg.pf_location(2) * (ny - pf_ny);
pf_kz = (1:pf_nz) + arg.pf_location(3) * (nz - pf_nz);

% ------- begin Xu algorithm -----------

% step 1: initial image recon based on zero fill
kspace_init = zeros(nx, ny, nz); % zero fill
if arg.fill_conj % fill other "half" of k-space with complex conjugate
		% (useful if image is real, but probably not otherwise)
	flip = @(i,n) mod(-(i-1), n) + 1;
	kspace_init(flip(pf_kx,nx), flip(pf_ky,ny), flip(pf_kz,nz)) = ...
		conj(partial_kspace);
end
kspace_init(pf_kx, pf_ky, pf_kz) = partial_kspace;
img_est = ifftn(ifftshift(kspace_init));

if 0 % check boundary of initial k-space
	im plc 1 2
	im_log(1, partial_kspace, 'init')
	im_log(2, kspace_init, 'init')
	keyboard
end

if ~isempty(arg.init)
	img_est = arg.init; % user-provided initial image
end

if arg.show
	im(3, abs(img_est)), title '|img_est|'
	im(8, angle(img_est), [-1 1]*pi), title '< img_est'
	im_log(13, kspace_init), title 'kspace init'
end

% step 2 and 3: Hanning tapered central k-space
lf_kspace = zeros(nx, ny, nz);
lf_kspace(lf_kx, lf_ky, lf_kz) = kspace_init(lf_kx, lf_ky, lf_kz);

% taper edges of central k-space
if arg.window_step3 > 0
	hanning_dims_step3 = [numel(lf_kx) numel(lf_ky) numel(lf_kz)];
	hanning_step3 = ir_mri_pf_gen_window( ...
		arg.window_step3, hanning_dims_step3, is_3d); % "filter"
	lf_kspace(lf_kx, lf_ky, lf_kz) = ...
		lf_kspace(lf_kx, lf_ky, lf_kz) .* hanning_step3;
end

% step 4: estimate phase from lf image
img_lf = ifftn(ifftshift(lf_kspace));
% img_lf = reale(img_lf); % use only for testing ideal real case
phi = angle(img_lf);

% build apodizer for step 8 outside loop
if arg.window_step8 > 0
	hanning_dims_step8 = [pf_nx pf_ny pf_nz] + arg.window_step8;
	hanning_step8 = ir_mri_pf_gen_window(...
		arg.window_step8, hanning_dims_step8, is_3d);
%	im(10, hanning_step8)

	if is_3d
		hanning_step8_crop = hanning_step8( ...
			(1:pf_nx) + arg.window_step8*(1-arg.pf_location(1)), ...
			(1:pf_ny) + arg.window_step8*(1-arg.pf_location(2)), ...
			(1:pf_nz) + arg.window_step8*(1-arg.pf_location(3)));
	else
		hanning_step8_crop = hanning_step8(...
			(1:pf_nx) + arg.window_step8*(1-arg.pf_location(1)), ...
			(1:pf_ny) + arg.window_step8*(1-arg.pf_location(2)));
	end
%	im(15, hanning_step8_crop)

	hanning_step8_full = zeros(nx, ny, nz);
	hanning_step8_full(pf_kx, pf_ky, pf_kz) = hanning_step8_crop;
else
%	hanning_step8_full = abs(kspace_init) > 0; % no - dangerous!
	hanning_step8_full = ones(size(kspace_init));
end

if arg.show % show everything before loop begins
	ktick = @(k, n) unique([1 n minmax(k)']);
	im(12, hanning_step8_full), title 'hanning step8'
	xtick(ktick(pf_kx, nx))
	ytick(ktick(pf_ky, ny))
	im(4, abs(img_lf)), title '|img_lf|'
	im(9, phi, [-1 1]*pi), title '< img_lf'
	im_log(14, lf_kspace), title 'lf_kspace'
	xtick(ktick(lf_kx, nx))
	ytick(ktick(lf_ky, ny))
	im(11, hanning_step3), title 'hanning_step3'
end

% pocs: image-domain phase constraint and k-space data constraint
for iter = 1:arg.niter
	% step 5:
	rho1 = abs(img_est) .* exp(1i*phi);
	% note: if img_est is real but has negative values, this abs
	% will eliminate those negative values, changing the spectrum!

	if arg.chat
	% examine how much it has changed
		change2 = norm(col(rho1-img_est)) / norm(col(img_est));
		changei = max(abs(col(rho1-img_est))) / max(abs(col(img_est)));
		printm('change = %4.2f%% %4.2f%%', change2 * 100, changei * 100)
	end

	% step 6:
	s1 = fftshift(fftn(rho1));

	% step 7 & 8: enforce data consistency but smoothly
	full_kspace = hanning_step8_full .* kspace_init ...
		+ (1-hanning_step8_full) .* s1;

	% step 9:
	img_est = ifftn(ifftshift(full_kspace));

	if arg.show
	%	phi_diff = phi - angle(img_est);
	%	im_log(9, fftshift(fftn(exp(1i*phi_diff))))
		im(5, abs(img_est)), title '|img_est|'
		im(10, angle(img_est)), title '\angle < img_est'
	%	im_log(14, s1), title 's1 spectrum'
		im_log(15, full_kspace), title 'new spectrum'
		drawnow
%	prompt
	end
end

end % ir_mri_partial_fourier_3d


% ir_mri_pf_gen_window()
function hanning_window = ir_mri_pf_gen_window(window_length, dims, is_3d)
% for apodization, window length refers to transition band length
% dims is size of window
if any(dims(1:(end-(1-is_3d))) <= window_length*2)
	fail('transition band wider than Hanning window dims')
	keyboard
end

% add extra 3: 2 because end points are zero, 1 in the middle for expanding
hann_len = (window_length+1)*2 + 1; % always odd
if ir_is_octave % trick: octave 'hanning' = matlab 'hann'
	hanning_1D = col(hanning(hann_len));
else % matlab
	hanning_1D = col(hann(hann_len));
end
%minmax( hanning_1D ) % [0 1]
hanning_2D = hanning_1D * hanning_1D';
if is_3d
	hanning_3D = repmat(hanning_2D, [1 1 hann_len]) .* ...
	repmat(permute(hanning_1D, [2 3 1]), [hann_len hann_len 1]);
	hanning_crop = hanning_3D(2:end-1,2:end-1,2:end-1) ./ ...
		max(col(hanning_3D));
else
	hanning_crop = hanning_2D(2:end-1,2:end-1) ./ max(col(hanning_2D));
end

% hanning_stretch = cat(1, hanning_crop(1:window_length,:,:), repmat(hanning_crop(window_length+1,:,:), dims(1) - 2*window_length,1,1), hanning_crop(window_length+2:end,:,:));

tmp = hanning_crop(window_length+1,:,:);
tmp = repmat(tmp, [dims(1) - 2*window_length, 1, 1]);
hanning_stretch = cat(1, hanning_crop(1:window_length,:,:), ...
	tmp, hanning_crop(window_length+2:end,:,:));

% hanning_stretch_2D = cat(2, hanning_stretch(:,1:window_length,:), repmat(hanning_stretch(:,window_length+1,:),1, dims(2) - 2*window_length,1), hanning_stretch(:,window_length+2:end,:));

tmp = hanning_stretch(:, window_length+1, :);
tmp = repmat(tmp, [1, dims(2) - 2*window_length, 1]);
hanning_window = cat(2, hanning_stretch(:,1:window_length,:), ...
	tmp, hanning_stretch(:,window_length+2:end,:));

if is_3d
	tmp1 = hanning_window(:,:,1:window_length),
	tmp2 = hanning_window(:,:,window_length+1);
	tmp3 = repmat(tmp2, [1, 1, dims(3) - 2*window_length]);
	hanning_window = cat(3, tmp1, tmp3, ...
		hanning_window(:,:,window_length+2:end));
end

end % ir_mri_pf_gen_window()


function ir_mri_partial_fourier_3d_test
end % test
