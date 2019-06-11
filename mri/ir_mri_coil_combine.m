 function [xhat, smap, ref] = ir_mri_coil_combine(ykj, varargin)
%function [xhat, smap, ref] = ir_mri_coil_combine(ykj, varargin)
%|
%| Given multiple MRI surface coil images, estimate their sensitivity maps,
%| relative to the sqrt-sum-of-squares (SSoS) map, and then perform the
%| "best" linear combination using the sensitivity maps.
%| todo: currently ignores noise correlations
%|
%| This routine provides an alternative to the usual "sqrt-sum-of-squares"
%| approach to combining multiple coil data.  It produces a complex-valued
%| output image that has approximately a complex gaussian distribution at
%| each voxel.  However, it does not preserve image phase information
%| because the image phase is included in the sensitivity map phase estimates.
%| So it might be logical to just take the real part of the output image.
%| The output is "approximately" gaussian only because of random errors in
%| the sensitivity maps estimated within.
%|
%| The variance is space-dependent because of coil sensitivity variations.
%| todo: return variance map!
%|
%| Signal model: y_kj = s_kj f_j + noise_kj
%|	s_kj	sensitivity map (relative to SSoS reference image)
%|	f_j	unknown underlying object
%|	k: coil index
%|	j: voxel index
%|
%| in
%|	ykj	[(N) ncoil]	noisy complex images (2D or 3D) for each coil
%|
%| options
%|	bodycoil [(N)]		reference image (optional) (default: SSoS)
%|	thresh0			zero reference image values below this fraction
%|				 of peak value (default: 0.05)
%|	thresh1			fraction of reference peak used for median
%|				 initial value in "background" (default: 0.05)
%|
%| out
%|	xhat	[(N)]		coil combination image
%|	smap	[(N) ncoil]	estimated sensitivity maps
%|
%| Uses the 'chol' (Cholesky) option of mri_sensemap_denoise.m which is
%| fast for medium sized images but may use much memory for large images.
%|
%| Copyright 2014-08-18, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(ykj, 'test')
	ir_mri_coil_combine_test
return
end

% defaults
arg.thresh0 = 0.05; % below this set reference image to 0
arg.thresh1 = 0.05; % below this set reference image to median
arg.bodycoil = [];
arg.chol = true;

arg = vararg_pair(arg, varargin);

tmp = size(ykj);
NN = tmp(1:end-1);
ncoil = tmp(end);

% default reference image
if isempty(arg.bodycoil)
	ref = ir_mri_coil_combine_ref(ykj, arg.thresh0);
else
	ref = arg.bodycoil;
end

smap = mri_sensemap_denoise(ykj, 'thresh', arg.thresh1, ...
	'bodycoil', ref, 'chol', 1, 'niter', 1);

coildim = ndims(ykj);
numer = sum(conj(smap) .* ykj, coildim);
denom = sum(abs(smap).^2, coildim);
xhat = div0(numer, denom); % LS solution at each pixel


% ir_mri_coil_combine_ref()
function ref = ir_mri_coil_combine_ref(ykj, thresh0)
coildim = ndims(ykj);
ref = sqrt(sum(abs(ykj).^2, coildim)); % sqrt-sum-of-squares
peak = max(abs(ref(:)));
ref(abs(ref) < peak * thresh0) = 0; % zero low SNR pixels
%ref = ref / peak; % normalize to 1


% ir_mri_coil_combine_test()
% built-in test/example
function ir_mri_coil_combine_test

if 1 % xtrue
	f.dir = [path_find_dir('mri') '/../data/mri/'];
	f.xtrue = [f.dir 'brainweb_t1.jpg'];
	xtrue = single(imread(f.xtrue)');
	xtrue = xtrue(2:end-1,2:end-1); % make it 256^2
	xtrue = downsample2(xtrue, 4); % now 64^2
	[nx, ny] = size(xtrue);
	atrue = 2*pi * (-0.5+([1:nx]'/nx * [1:ny]/ny).^0.5); % smooth phase
	xtrue = xtrue .* exp(1i * atrue); % phase
%	im clf, im('hsv', angle(xtrue), [-pi pi]), cbar, return

	im plc 3 4
	im(1, abs(xtrue), '$|x|$'), cbar
	im(2, 'hsv', angle(xtrue), [-pi pi], '$\angle x$'), cbar
	mask = conv2(abs(xtrue), ones(5), 'same') > 0;
%	im(11, abs(xtrue) + 100*mask), return % check mask
end

if 1 % smap
	ncoil = 4;
	mask4 = repmat(mask, [1 1 ncoil]);
	smap = ir_mri_sensemap_sim('nx', nx, 'ny', ny, 'dx', 192/nx, ...
		'ncoil', ncoil, 'rcoil', 100);
%	smap = smap(:, :, [ncoil 1:(ncoil-1)]); % phase of last coil is cooler
	tmp = sqrt(sum(abs(smap).^2, 3));
	smap = smap / tmp(end/2,end/2); % normalize true smap for simplicity
	im(5, mask4 .* abs(smap), '$|s|$'), cbar
	im(6, 'hsv', mask4 .* angle(smap), [-pi pi], '$\angle s$'), cbar
	im(12, mask .* sqrt(sum(abs(smap).^2, 3)), 's SSoS'), cbar
end

if 1 % y
	ytrue = smap .* repmat(xtrue, [1 1 4]);
	rng(0)
	snr2sigma = @(db, yb) exp(-db/20) * norm(yb(:)) / sqrt(numel(yb)) / sqrt(2); % for complex noise
	sig = snr2sigma(50, ytrue);
%	sig = 0; % noiseless
	ykj = ytrue + sig * (randn(size(ytrue)) + 1i * randn(size(ytrue)));
	im(7, abs(ykj), '$|y|$'), cbar
	im(8, 'hsv', ir_unwrap(angle(ykj)), [-pi pi], '$\angle y$'), cbar
end

	[xhat shat ref] = ir_mri_coil_combine(ykj);
	ssos = sqrt(sum(abs(ykj).^2, 3));
	xhat = xhat * ir_best_scale(xhat, abs(xtrue));
	ssos = ssos * ir_best_scale(ssos, abs(xtrue));

	pr nrms(xhat(:), col(abs(xtrue)))
	pr nrms(ssos(:), col(abs(xtrue)))
%	pr corrcoef(abs(xtrue(:)), ssos(:))
%	pr corrcoef(abs(xtrue(:)), real(xhat(:)))

%	xhat = reale(xhat);
%	im(3, abs(xhat), '|x hat|'), cbar
%	im(4, angle(xhat), [-pi pi], '< x hat'), cbar
	im(3, real(xhat), 'Re($\hat{x}$)'), cbar
	im(4, imag(xhat), 'Im($\hat{x}$)'), cbar
	im(9, mask4 .* abs(shat), '$|\hat{s}|$'), cbar
	im(10, 'hsv', mask4 .* angle(shat), [-pi pi], '$\angle \hat{s}$'), cbar
	im(11, ssos, 'x SSoS'), cbar
