% sense_homodyne_recon_demo.m
% demos use of sense_homodyne_recon.m, which combines
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
sense_direction = 2;
homodyne_direction = 2;

% turn figures off for SENSE function
figs = 0;

% number of coils
nc = 4; 
% factor of undersampling
np = 2; 

% amount of overlap for homodne reconstruction
overlap = 10;
regularizer = 'none'; % can also be 'tikhonov'

% read brain image
fsep = filesep;
img_dir = [path_find_dir('mri') fsep '..' fsep 'data' fsep 'mri'];
img_path = [img_dir fsep 'brainweb_t1.jpg'];
image = double(imread(img_path));
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
im_fft = [];
reduced_ffts = [];
sense_only_ffts = [];
for ii = 1:nc
    im_fft(:,:,ii) = fft2((mapped_im(:,:,ii)));
    if (sense_direction == 1)
        reduced_fft = im_fft(1:np:end,:,ii);
    else
        reduced_fft = im_fft(:,1:np:end,ii);
    end
    sense_only_ffts(:,:,ii) = reduced_fft;
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
sense_only_ffts = sense_only_ffts + sigma*randn(size(sense_only_ffts));

%% construct mask
% only have access to reduced_ffts (after SENSE and homodyne reductions)
% aliased_ims = (ifft2(ifftshift(reduced_ffts)));
% figure; imagesc(real(aliased_ims(:,:,1)));
mask = any(abs(mapped_im)>0.005,3);

%% 
[demod_recon_im, recon_im, lp_im] = sense_homodyne_recon(reduced_ffts, smap, overlap, homodyne_direction, sense_direction, regularizer, nc, np, red_dims, mask, figs);
sense_only_recon = sense_recon(smap, sense_only_ffts, sense_direction, np, regularizer, mask, figs);

%% plot SENSE & homodyne reconstructed images

figure; subplot(331); imshow(abs(image),[]); colorbar; 
title('original image');

subplot(332); imshow(abs(sense_only_recon),[]); colorbar;
title('abs(SENSE reconstructed image)');

subplot(333); imshow(abs(image)-abs(sense_only_recon),[]); colorbar;
sense_only_MSE = sum(sum((abs(image)-abs(sense_only_recon)).^2))/prod(dims)
title(['diff b/n SENSE only and orig, MSE: ' num2str(sense_only_MSE)]);

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

