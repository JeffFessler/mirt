 function [wmap, wconv] = mri_field_map_reg(yik, etime, varargin)
%function [wmap, wconv] = mri_field_map_reg(yik, etime, varargin)
%|
%| MRI field map estimation using regularization and optimization transfer,
%| based on a multiple image model:
%|	y_ik = mag_i * phase_ik * decay_ik + noise_ik
%|	phase_ik = exp(1i * wmap_i * (etime_k - etime_1))
%|	decay_ik = exp(-R2_i * (etime_k - etime_1))
%|	mag_i	unknown magnitude map
%|	wmap_i	unknown field map
%|	R2_i	unknown relaxation map
%| and etime_k are echo times (or relative echo times). Often etime(1)=0.
%| Units of wmap will be reciprocal of units of etime_k vector.
%| So if etime has units seconds, then wmap has units rad/sec.
%|
%| in
%|	yik	[(N) nset]	complex images at different echo times
%|				(these should be undistorted)
%|				(N) can be (nx,ny) or (nx,ny,nz) etc.
%|	etime	[nset]		vector of echo times [usually in seconds]
%|
%| options
%|	l2b		log_2(beta), regularization parameter (default: -3)
%|	order		regularization order (default: 2)
%|	niter		# of iterations (default: 40)
%|	fmax		percent of median used in calculating scale factor (.1)
%| 	wthresh		percent of median used in masking out for RMSerror (.1)
%|	winit	[(N)]		initial field map for iterating
%|				(default: estimate from first two scans)
%|	mask	[(N)]		logical support mask (only estimate within this)
%|
%| out
%|	wmap	[(N)]	regularized estimate of field map [rad/sec]
%|	wconv	""	conventional estimate from first two scans
%|
%| The first data set (usually) has no phase offset (etime(1) = 0).
%| The second data set (usually) tries to have no phase wrapping,
%| i.e., etime(2)-etime(1) usually is small.
%| The other (optional) datasets (usually) have phase wrapping, but higher SNR.
%|
%| This algorithm is based on the paper:
%| A K Funai, J A Fessler, D T B Yeo, V T Olafsson, D C Noll
%| "Regularized field map estimation in MRI"
%| IEEE Trans. Med. Imaging, 27(10):1484-94, Oct. 2008.
%| http://dx.doi.org/10.1109/TMI.2008.923956
%|
%| See the built-in test routine for a self-contained example.
%|
%| Caution: the sign of the field map estimated here is the opposite (negative)
%| of the sign of the field map needed for input to the Gmri object.
%|
%| Copyright 2007-12-15, Amanda Funai & Jeff Fessler, University of Michigan

if nargin >= 1 && streq(yik, 'test')
	wmap = mri_field_map_reg_test(yik, varargin{:});
	if ~nargout, clear wmap, end
return
end

if nargin < 2, help(mfilename), error(mfilename), end

% Specify defaults & check argument inputs
arg.l2b = -3;
arg.order = 2;
arg.niter = 40;
arg.niter_init = 10;
arg.fmax = 0.1;
arg.wthresh = 0.1;
arg.winit = [];
arg.mask = [];

arg = vararg_pair(arg, varargin);

% find (N) aka dim_im
dim_yik = size(yik);
if length(dim_yik) == 2
	dim_im = [dim_yik(1) 1];
else
	dim_im = dim_yik(1:(end-1));
end
nset = dim_yik(end);

if size(etime,2) ~= nset
	fail('etime vector must be length %d', nset)
end

if isempty(arg.mask)
	arg.mask = true(dim_im);
else
	jf_equal(size(arg.mask), dim_im)
end

%% scale/normalize data to make regularization effect invariant to scaling factors
% todo: more testing needed here
if 1
	[yik yi2] = mri_field_map_reg_scale(yik, etime, arg.fmax);
else
	yik = ir_mri_field_map_reg_scale(yik, etime, 'fmax', arg.fmax);
	yi2 = ir_mri_field_map_reg_scale(stackpick(yik,1:2), ...
		etime(1:2), 'fmax', arg.fmax);
%	im plc 2 2, im(1, yik), im(2, yi2), im(3, tmp1), keyboard
end

% columnize/mask y for use in SPS and in calculating cost
ycolk = reshape(yik, prod(dim_im), nset);	ycolk = ycolk(arg.mask,:);
ycol2 = reshape(yi2, prod(dim_im), 2);		ycol2 = ycol2(arg.mask,:);

% iterative estimation with quadratic regularization
Rq = Reg1(arg.mask, 'beta', 2^arg.l2b, 'order', arg.order);

% if needed, use regularized initialization from first two data sets
if isempty(arg.winit)
	% start with conventional estimate
	wconv = angle(stackpick(yi2,2) .* conj(stackpick(yi2,1))) ...
		/ (etime(2) - etime(1));
	arg.winit = unwrapping_sps_manysets(wconv(arg.mask), ycol2, ...
		etime(1:2), Rq, 'niter', arg.niter_init);
	arg.winit = embed(arg.winit(:,end), arg.mask);
%	im plc 1 2, im(1, arg.winit), im(2, isnan(arg.winit)), prompt
	mag1 = abs(stackpick(yik,1));
	good = mag1 > arg.wthresh * max(mag1(:));
	arg.winit(~good) = mean(arg.winit(good));
else
	jf_equal(size(arg.mask), arg.winit)
end

if any(isnan(arg.winit(:)))
	im clf, im(arg.winit)
	warn 'bug: nan in arg.winit'
	keyboard
end

% estimate field map from all scans using the given initialization
wmap = unwrapping_sps_manysets(arg.winit(arg.mask), ycolk, etime, Rq, ...
	'niter', arg.niter);
wmap = embed(wmap(:,end), arg.mask);


% mri_field_map_reg_scale()
% in
%	yik	[(N) nset]
%	etime	[nset]
% out
%	yik	[(N) nset]
%	yi2	[(N) 2]		(first two echos only)
%
function [yik yi2] = mri_field_map_reg_scale(yik, etime, fmax)

% todo: replace with ir_mri_field_map_reg_scale after checking compatability

dim_yik = size(yik);
yik = reshapee(yik, [], dim_yik(end)); % [*N nset]

[nn nset] = size(yik);

% Scale by median of first set of data to get rid of large mag_j effects
y1 = abs(yik(:,1));
scalefactor = median(y1(y1(:) > fmax * max(y1(:))));
if scalefactor == 0
	fail 'median is zero?'
end
yik = yik / scalefactor;

% Try to compensate for R2 effects on effective regularization.
% (This can be slightly different scaling for each pixel.)
% "approximate oracle scaling"

wjfactor = zeros(nn,1);
wjfactor2 = zeros(nn,1);
for j=1:nset
	for k=1:nset
		tmp = yik(:,j) .* yik(:,k) * (etime(k)-etime(j));
		wjfactor = wjfactor + abs(tmp).^2;
		if j <= 2 && k <= 2
			wjfactor2 = wjfactor2 + abs(tmp).^2;
		end
	end
end
totalsum = sum(abs(yik).^2,2) .* abs(yik(:,1)).^2;
wjfactor = sqrt(div0(wjfactor, totalsum));

totalsum2 = div0(abs(yik(:,2)).^2 + abs(yik(:,1)).^2, abs(yik(:,1)).^2);
wjfactor2 = sqrt(div0(wjfactor2, totalsum2));

% Scale data by this factor
yi2 = [div0(yik(:,1), wjfactor2), div0(yik(:,2), wjfactor2)];
for j=1:nset
	yik(:,j) = div0(yik(:,j), wjfactor);
end

yik = reshape(yik, dim_yik); % [(N) nset]
yi2 = reshape(yi2, [dim_yik(1:end-1) 2]); % [(N) 2]


% mri_field_map_reg_test()
% built-in test/example
function wmap = mri_field_map_reg_test(type, varargin)

% simulate multiple images based on a given mag & field map.
printm 'simulate noisy data'
etime = [0 2 10] * 1e-3; % echo times
nset = length(etime);
R2 = 20; % 1/sec decay
SNR = 20; % dB

dir = [path_find_dir('mri') '/phase-data/'];
wtrue = 2*pi * fld_read([dir 'fieldmap128.fld']); % "true" fieldmap
mag = fld_read([dir 'mag128.fld']); % "true" magnitude
[nx ny] = size(mag); % 128

im plc 2 4
im(1, wtrue/(2*pi), 'true field map', [-40, 128]), cbar('Hz')
im(5, mag, 'true mag'), cbar

image_power = 10*log10(sum(sum(mag.^2))/(nx*ny));
noise_power = image_power - SNR;
noise_std = sqrt(10^(noise_power/10));
noise_std = noise_std / 2; % because complex

yik = zeros(nx,ny,nset);
for kk=1:nset
	yik(:,:,kk) = mag ...
		.* exp(1i * wtrue * (etime(kk) - etime(1))) ...
		.* exp(-R2 * (etime(kk) - etime(1)));
end
rng(0)
yik = yik + noise_std * (randn(size(yik)) + 1i * randn(size(yik)));

%yik(1:10,1:10,:) = 0; % stress test all 0 regions

printm 'estimate field map'
[wmap wconv] = mri_field_map_reg(yik, etime, 'l2b', -6, varargin{:});

mask = mag > 0.05 * max(mag(:));
im(8, mask), cbar
titlef('Mask for RMSE')

im(2, wconv / (2*pi), 'Conventional field map', [-40,128]), cbar('Hz')
err = (wconv - wmap) / (2*pi);
im(6, err, 'error', [-25 25]), cbar('Hz')
xlabelf('RMSE = %.1f Hz', sqrt(mean(err(mask).^2)))

im(3, wmap / (2*pi), 'Regularized field map', [-40,128]), cbar('Hz')
err = (wtrue - wmap) / (2*pi);
im(7, err, 'error', [-25 25]), cbar('Hz')
xlabelf('RMSE = %.1f Hz', sqrt(mean(err(mask).^2)))
