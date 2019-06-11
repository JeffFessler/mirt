% SENSE_homodyne_recon_demo.m
% demos use of SENSE_homodyne_recon.m, which combines
% SENSE and partial Fourier reconstruction. Input image can be rectangular,
% but both dimensions should be even
%
% resulting figure shows the improvement from demodulating with the
% reference image and compares MSE to reconstruction from an image that
% only undergoes SENSE undersampling and not homodyne
%
% 2012-06-18 Mai Le

% can undersample in either direction for either method
% values should be {1,2}
SENSE_direction = 2;
homodyne_direction = 1;

% turn figures off for SENSE function
figs = 0;

% number of coils
nc = 6; 
% factor of undersampling
np = 4; 

% amount of overlap for homodne reconstruction
overlap = 10;
regularizer = 'none'; % can also be 'tikhonov'

% read brain image
image = double(imread('./data/mri/brainweb_t1.jpg','jpg'));
% requires image to be multiple of np in dimension of undersampling
if (SENSE_direction == 1)
    image = image(1:floor(size(image,1)/np)*np,:);
else
    image = image(:,1:floor(size(image,2)/np)*np);
end
dims = size(image);

% figure;
smap = mri_sensemap_sim('nx',dims(1),'ny',dims(2),'ncoil',nc);

% simulate linear phase 
[xx,yy] = ndgrid(-dims(1)/2:dims(1)/2-1,-dims(2)/2:dims(2)/2-1);
ph_max = pi/5;
ph = ph_max*xx/(dims(1)/2)+ph_max*yy/(dims(2)/2); 
image = image.*exp(1i*ph);

% replicate image
im_rep = repmat(image,[1 1 nc]);

% apply sensitivity maps
mapped_im = im_rep.*smap;

%% Generate partial k-space data
im_fft = zeros(dims(1),dims(2),nc);
reduced_ffts = [];
SENSE_only_ffts = [];
for ii = 1:nc
    im_fft(:,:,ii) = fft2((mapped_im(:,:,ii)));
    if (SENSE_direction == 1)
        reduced_fft = im_fft(1:np:end,:,ii);
    else
        reduced_fft = im_fft(:,1:np:end,ii);
    end
    SENSE_only_ffts(:,:,ii) = reduced_fft;
    reduced_fft = fftshift(reduced_fft); 
    red_dims = size(reduced_fft);
    if (homodyne_direction == 1)
        reduced_fft = reduced_fft(1:min(red_dims(1),ceil(red_dims(1)/2)+overlap+1),:);
    else 
        reduced_fft = reduced_fft(:,1:min(red_dims(2),ceil(red_dims(2)/2)+overlap+1));
    end
    reduced_ffts(:,:,ii) = reduced_fft;
end

%% % introduce Gaussian noise
sigma = 0;
reduced_ffts = reduced_ffts + sigma*randn(size(reduced_ffts));
SENSE_only_ffts = SENSE_only_ffts + sigma*randn(size(SENSE_only_ffts));

%% construct mask
mask = any(abs(mapped_im)>0.005,3);

%% 
[demod_recon_im, recon_im, lp_im] = SENSE_homodyne_recon(reduced_ffts, smap, overlap, homodyne_direction, SENSE_direction, regularizer, nc, np, red_dims, mask, figs);
SENSE_only_recon = SENSE_recon(smap, SENSE_only_ffts, SENSE_direction, np, regularizer, mask, figs);

%% plot SENSE & homodyne reconstructed images

figure; subplot(331); imshow(abs(image),[]); colorbar; 
title('original image');

subplot(332); imshow(abs(SENSE_only_recon),[]); colorbar;
title('abs(SENSE reconstructed image)');

subplot(333); imshow(abs(image)-abs(SENSE_only_recon),[]); colorbar;
SENSE_only_MSE = sum(sum((abs(image)-abs(SENSE_only_recon)).^2))/prod(dims)
title(['diff b/n SENSE only and orig, MSE: ' num2str(SENSE_only_MSE)]);

subplot(334); imshow(angle(recon_im),[]); colorbar;
title('phase of recon w/o demodulation');

subplot(335); imshow(real(recon_im),[]); colorbar;
title('real of recon w/o demodulation');

subplot(336); imshow(abs(image)-real(recon_im),[]); colorbar;
nondemod_MSE = sum(sum((abs(image)-real(recon_im)).^2))/prod(dims)
title(['diff b/n nondemod SENSE&h and orig, MSE: ' num2str(nondemod_MSE)]);

subplot(337); imshow(angle(lp_im),[]); colorbar;
title('phase of reference image');

subplot(338); imshow((demod_recon_im),[]); colorbar;
title('SENSE & homodyne, demod reconstructed image');

subplot(339); imshow(abs(image)-demod_recon_im,[]); colorbar;
demod_MSE = sum(sum((abs(image)-demod_recon_im).^2))/prod(dims)
title(['diff b/n demod SENSE&h and orig, MSE: ' num2str(demod_MSE)]);

