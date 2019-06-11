  function [fwhm_best, costs, im_best] = ...
	fwhm_match(true_image, blurred_image, fwhms)
%|function [fwhm_best, costs, im_best] = ...
%|	fwhm_match(true_image, blurred_image, fwhms)
%|
%| given a blurred_image of a true_image, find the FHWM of a Gaussian kernel
%| that, when convolved to the true_image, yields the smoothed image
%| that best matches blurred_image.
%|
%| the set of FWHM values given in the array fwhms is tried.
%|
%| Copyright 2001-8-30, Jeff Fessler, University of Michigan

if nargin == 1 && streq(true_image, 'test'), fwhm_match_test, return, end
if nargin < 2, ir_usage, end

if nargin < 3
	fwhms = 0:0.5:4;
end

costs = zeros(size(fwhms));
cost_min = Inf;
for ii=1:length(fwhms)
	fwhm = fwhms(ii);
	kern = gaussian_kernel(fwhm);
	psf = kern * kern';
	tmp = conv2(true_image, psf, 'same');
	costs(ii) = norm(tmp(:) - blurred_image(:)) / norm(true_image(:));
	if costs(ii) < cost_min
		im_best = tmp;
	end
end

[dummy ibest] = min(costs);
if ibest == 1 || ibest == length(fwhms)
	warning 'need wider range of fwhms'
end
fwhm_best = fwhms(ibest);


% fwhm_match_test
function fwhm_match_test

% pyramidal PSF to stress the approach
psf1 = [0:5 4:-1:0]; psf1 = psf1 / sum(psf1); psf = psf1' * psf1;
true_image = zeros(128); true_image(64:96,64:96) = 1;
blurred_image = conv2(true_image, psf, 'same');
im plc 2 2
im(1, true_image, 'True Image')
im(2, blurred_image, 'Blurred Image')

fwhms = [2:0.25:8];
[fwhm_best, costs] = fwhm_match(true_image, blurred_image, fwhms);
np = length(psf);	ip = -(np-1)/2:(np-1)/2;
kern = gaussian_kernel(fwhm_best);
nk = length(kern);	ik = -(nk-1)/2:(nk-1)/2;

if im
	im subplot 3
	plot(fwhms, costs, 'c-o', fwhm_best, min(costs), 'yx')
	xlabel FWHM, ylabel Cost, title 'Cost vs FWHM' 
	im subplot 4
	plot(ip, psf1, '-o', ik, kern(:), '-+')
	xlabel pixel, title 'PSF profile: actual and Gaussian fit'
end
