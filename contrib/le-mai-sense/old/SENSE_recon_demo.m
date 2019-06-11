% SENSE_recon_demo.m
% - demonstrates use of SENSE_recon.m, which implements the SENSE
% reconstruction
% - figure shows how Tikhonov regularization mitigates the increase in error
% as noise variance increases 
% - applies complex AWGN to kspace
% 
% 2012-06-08, Mai Le


% number of coils
nc = 4; 
% factor of undersampling
np = 2;
% dimension to undersample
reduced_dim = 1;
regularizer1 = 'none';
regularizer2 = 'tikhonov';
figs = 0;

image = double(imread('./data/mri/brainweb_t1.jpg','jpg'));

dims = size(image);
% introduce planar phase
[xx,yy] = ndgrid(-dims(1)/2:dims(1)/2-1,-dims(2)/2:dims(2)/2-1);
ph_max = pi/5;
ph = ph_max*xx/(dims(1)/2)+ph_max*yy/(dims(2)/2);
image = image.*exp(1i*ph);

% generate sensitivity maps
smap = mri_sensemap_sim('nx',dims(1),'ny',dims(2),'ncoil',nc);

im_rep = repmat(image,[1 1 nc]);

% apply sensitivity maps
mapped_im = im_rep.*smap;

% throw away k-space data
im_fft = [];
reduced_fft = [];
for ii = 1:nc
    im_fft(:,:,ii) = fft2(mapped_im(:,:,ii));
    if (reduced_dim == 1)
        reduced_fft(:,:,ii) = im_fft(1:np:end,:,ii);
    else
        reduced_fft(:,:,ii) = im_fft(:,1:np:end,ii);
    end
end
%% perform SENSE reconstruction for varying levels of noise
sigmas = 0:25:150;
for ii = 1:length(sigmas)
    % introduce complex AWGN
    sigma = sigmas(ii);
    noisy_reduced_fft = reduced_fft + sigma*randn(size(reduced_fft)) + i*sigma*randn(size(reduced_fft));
    
    % SENSE reconstruction
    recon_im1(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, reduced_dim, np, regularizer1, figs);
    recon_im2(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, reduced_dim, np, regularizer2, figs);
    
    error1(ii) = sqrt((sum(sum((real(image)-real(recon_im1(:,:,ii))).^2))+sum(sum((imag(image)-imag(recon_im1(:,:,ii))).^2)))/prod(dims));
    error2(ii) = sqrt((sum(sum((real(image)-real(recon_im2(:,:,ii))).^2))+sum(sum((imag(image)-imag(recon_im2(:,:,ii))).^2)))/prod(dims)) ;
end
%% plot error as a function of noise variance
figure; plot(sigmas.^2, error1,'k');
hold on; plot(sigmas.^2, error2, 'r'); legend('no regularization', 'tikhonov');
xlabel('variance of complex AWGN in kspace');
ylabel('RMS error in reconstructed images');
title('effect of Tikhonov regularization on noisy SENSE reconstruction')
%% plot reconstructed image at sigma(ii)
ii = 4; % change as desired
figure; subplot(331); imshow(real(image),[]); colorbar; 
title('real of original image');
subplot(332); imshow(imag(image),[]); colorbar; 
title(['imag of original image, noise variance: ' num2str(sigmas(ii).^2)]);
subplot(334); imshow(real(recon_im1(:,:,ii)),[]); colorbar;
title('real of recon image, no reg');
subplot(335); imshow(imag(recon_im1(:,:,ii)),[]); colorbar;
title('imag of recon image, no reg');
subplot(336); imshow(abs(real(image)-real(recon_im1(:,:,ii)))+abs(imag(image)-imag(recon_im1(:,:,ii))),[]); colorbar;
title('abs(differences), no reg');
subplot(337); imshow(real(recon_im2(:,:,ii)),[]); colorbar;
title('real of recon image, Tikhonov');
subplot(338); imshow(imag(recon_im2(:,:,ii)),[]); colorbar;
title('imag of recon image, Tikhonov');
subplot(339); imshow(abs(real(image)-real(recon_im2(:,:,ii)))+abs(imag(image)-imag(recon_im2(:,:,ii))),[]); colorbar;
title('abs(differences), Tikhonov');