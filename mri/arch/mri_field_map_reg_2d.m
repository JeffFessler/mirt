 function [wmap, wconv] = mri_field_map_reg_2d(yi, etime, varargin)
%function [wmap, wconv] = mri_field_map_reg_2d(yik, etime, varargin)
%
% MRI field map estimation using regularization and optimization transfer,
% based on a multiple image model:
%	y_ik = mag_i * phase_ik * decay_ik + noise_ik
%	phase_ik = exp(1i * wmap_i * (etime_k - etime_1))
%	decay_ik = exp(-R2_i * (etime_k - etime_1))
%	mag_i	unknown magnitude map
%	wmap_i	unknown field map
%	R2_i	unknown relaxation map
% and etime_k are echo times (or relative echo times).  Often etime(1)=0.
% Units of wmap will be reciprocal of units of etime_k vector.
% So if etime has units seconds, then wmap has units rad/sec.
%
% in
%	yik	[nx,ny,nsets]	complex images at different echo times
%				(these should be undistorted)
%	etime	[nsets]		vector of echo times
%
% options
%	l2b		log_2(beta), regularization parameter
%	order		regularization order (default: 2)
%	niter		# of iterations
%	fmax		percent of median used in calculating scale factor
% 	thresh		percent of median used in masking out for RMSerror
%	winit	[nx,ny]		initial field map for iterating
%				(default is to estimate from first two scans)
%	mask	[nx,ny]		logical support mask (only estimate within this)
%
% out
%	wmap	[nx,ny]	regularized estimate of field map	
%	wconv	""	conventional estimate from first two scans
%
% The first data set (usually) has no phase offset (etime(1) = 0).
% The second data set (usually) tries to have no phase wrapping,
% i.e., etime(2)-etime(1) usually is small.
% The other (optional) datasets (usually) have phase wrapping, but higher SNR.
%
% Copyright 2007-12-15, Amanda Funai & Jeff Fessler, The University of Michigan

if nargin >= 1 && streq(yi, 'test')
	wmap = mri_field_map_reg_2d_test(yi, varargin{:});
	if ~nargout, clear wmap, end
return
end

if nargin < 2, help(mfilename), error(mfilename), end

% Specify defaults & check argument inputs
arg.l2b = -3;
arg.order = 2;
arg.niter = 40;
arg.niter_init = 10;
arg.fmax = .1;
arg.wthresh = .1;
arg.winit = [];
arg.mask = [];

arg = vararg_pair(arg, varargin);

[nx,ny,nsets] = size(yi);

if size(etime,2) ~= nsets
	fail 'etime vector must be length nsets'
end

if isempty(arg.mask)
	arg.mask = true(nx,ny);
else
	jf_equal(size(arg.mask), [nx ny])
end

% scale data to make regularization effect invariant to scaling factors
[yi y2] = mri_field_map_reg_2d_scale(yi, etime, arg);

% Reshape y for use in SPS & in calculating cost
ycol = reshape(yi, nx*ny, nsets); ycol = ycol(arg.mask, :);
ycol2 = reshape(y2, nx*ny, 2); ycol2 = ycol2(arg.mask, :);

%% Iterative Estimation

Rq = Robject(arg.mask, 'type_denom', 'matlab', ...
	'beta', 2^arg.l2b, 'order', arg.order); % quadratic regularization

% If needed, use regularized initialization from first two data sets
if isempty(arg.winit)
	% start with conventional estimate
	wconv = angle(y2(:,:,2) .* conj(y2(:,:,1))) / (etime(2) - etime(1));
	arg.winit = unwrapping_sps_manysets(wconv(arg.mask), ycol2, ...
		etime(1:2), Rq, arg.niter_init);
	arg.winit = embed(arg.winit(:,end), arg.mask);
	mag1 = abs(yi(:,:,1));
	good = mag1 > arg.wthresh * max(mag1(:));
	arg.winit(~good) = mean(arg.winit(good));
end

% estimate field map from all scans using the given initialization
wmap = unwrapping_sps_manysets(arg.winit(arg.mask), ycol, etime, Rq, arg.niter);
wmap = embed(wmap(:,end), arg.mask);


%
% mri_field_map_reg_2d_scale()
%
function [yi y2] = mri_field_map_reg_2d_scale(yi, etime, arg);

[nx,ny,nsets] = size(yi);
y1 = yi(:,:,1);

% Scale by median of first set of data to get rid of large mag_j effects
scalefactor = median(y1(y1(:) > arg.fmax * max(y1(:))));
yi = yi / scalefactor;

% Try to compensate for R2 effects on effective regularization.
% (This can be slightly different scaling for each pixel.)
% "approximate oracle scaling"

wjfactor = zeros(nx,ny);
wjfactor2 = zeros(nx,ny);
for j=1:nsets
	for k=1:nsets
		tmp = yi(:,:,j) .* yi(:,:,k) * (etime(k)-etime(j));
		wjfactor = wjfactor + abs(tmp).^2;
		if j <= 2 && k <= 2
			wjfactor2 = wjfactor2 + abs(tmp).^2;
		end
	end
end
totalsum = sum(abs(yi).^2,3) .* abs(yi(:,:,1)).^2;
wjfactor = sqrt(wjfactor ./ totalsum);

totalsum2 = (abs(yi(:,:,2)).^2 + abs(yi(:,:,1)).^2) ./ (abs(yi(:,:,1)).^2);
wjfactor2 = sqrt(wjfactor2 ./ totalsum2);

% Scale data by this factor
y2(:,:,1) = yi(:,:,1) ./ wjfactor2;
y2(:,:,2) = yi(:,:,2) ./ wjfactor2;
for j=1:nsets
	yi(:,:,j) = yi(:,:,j) ./ wjfactor;
end


%
% built-in test/example
%
function wmap = mri_field_map_reg_2d_test(type, varargin)

% first make multiple datasets based on etime and a given mag & field map.
printm 'simulate noisy data'
etime = [0 2 10] * 1e-3;
nsets = length(etime);
R2 = 20; % 1/sec decay
SNR = 20; % dB

dir = [path_find_dir('mri') '/phase-data/'];
wtrue = 2*pi * fld_read([dir 'fieldmap128.fld']);
mag = fld_read([dir 'mag128.fld']);
[nx ny] = size(mag); % 128

im pl 2 3
im(1, wtrue/(2*pi), 'true field map', [-40, 128]), cbar('Hz')
im(4, mag, 'mag'), cbar

image_power = 10*log10(sum(sum(mag.^2))/(nx*ny));
noise_power = image_power - SNR;
noise_std = sqrt(10^(noise_power/10));
noise_std = noise_std / 2; % because complex 

yik = zeros(nx,ny,nsets);
for kk=1:nsets
	yik(:,:,kk) = mag ...
		.* exp(1i * wtrue * (etime(kk) - etime(1))) ...
		.* exp(-R2 * (etime(kk) - etime(1)));
end
rng(0)
yik = yik + noise_std * (randn(size(yik)) + 1i * randn(size(yik)));

printm 'estimate field map'
[wmap wconv] = mri_field_map_reg_2d(yik, etime, 'l2b', -6, varargin{:});

mask = mag > 0.05 * max(mag(:));
im(2, wconv / (2*pi), 'conventional field map', [-40,128]), cbar('Hz')
error = (wconv - wmap) / (2*pi);
im(5, error, 'error', [-25 25]), cbar('Hz')
xlabelf('RMS = %.1f Hz', sqrt(mean(error(mask).^2)))

im(3, wmap / (2*pi), 'regularized field map', [-40,128]), cbar('Hz')
error = (wtrue - wmap) / (2*pi);
im(6, error, 'error', [-25 25]), cbar('Hz')
xlabelf('RMS = %.1f Hz', sqrt(mean(error(mask).^2)))
