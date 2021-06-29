function [wmap, wconv] = mri_field_map_qm(yik, etime, varargin)
% This function is a wrapper for fmap_est_qm.m to make it callable in a
% similar fashion as mri_field_map_reg.m

%|	yik	[(N) nset ncoil]	complex images at different echo times
%|				(these should be undistorted)
%|				(N) can be (nx,ny) or (nx,ny,nz) etc.
%|	etime	[nset]		vector of echo times [usually in seconds]
%|
%| options
%|	see fmap_est_qm.m for options
%|
%| out
%|	wmap	[(N)]	regularized estimate of field map [rad/sec]
%|	wconv	""	conventional estimate from first two scans
%|
%| Written 2021-06-03, Melissa Haskell, University of Michigan

if nargin >= 1 && streq(yik, 'test')
	wmap = mri_field_map_qm_test(yik, varargin{:});
	if ~nargout, clear wmap, end
return
end

if nargin < 2, ir_usage, end

%% call function
func_name = 'fmap_est_qm';
[wmap, wconv] = call_fmap(yik, etime, func_name, varargin{:});



end



% mri_field_map_qm_test()
% built-in test/example, copied from mri_field_map_reg_test
function wmap = mri_field_map_qm_test(~, varargin)

% simulate multiple images based on a given mag & field map.
printm 'simulate noisy data'
etime = [0 2 10] * 1e-3; % echo times
nset = length(etime);
R2 = 20; % 1/sec decay
SNR = 20; % dB

wtrue = 2*pi * ir_get_data(fullfile('mri','2001-phase-data','fieldmap128.fld')); % "true" fieldmap
mag = ir_get_data(fullfile('mri','2001-phase-data','mag128.fld')); % "true" magnitude
[nx, ny] = size(mag); % 128

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
[wmap wconv] = mri_field_map_qm(yik, etime, 'l2b', -6, varargin{:});

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
end
