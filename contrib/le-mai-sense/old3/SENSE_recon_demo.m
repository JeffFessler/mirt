% SENSE_recon_demo.m
% - demonstrates use of SENSE_recon.m, which implements the SENSE
% reconstruction
% - figure shows how Tikhonov regularization mitigates the increase in error
% as noise variance increases 
% - applies complex AWGN to kspace
% - image must be even both directions
% 
% 2012-06-08, Mai Le

% comparison = 1 will demonstrate effect of masks and regularization on MSE,
% but demo takes much longer. comparison = 0 runs masked SENSE w/o
% regularization
comparison = 1;

% number of coils
nc = 4; 
% factor of undersampling
np = 3;
% dimension to undersample
reduced_dim = 2;
regularizer1 = 'none';
regularizer2 = 'tikhonov';
figs = 0;

image = double(imread('./data/mri/brainweb_t1.jpg','jpg'));
% requires image to be multiple of np in dimension of undersampling
if (reduced_dim == 1)
    image = image(1:floor(size(image,1)/np)*np,:);
else
    image = image(:,1:floor(size(image,2)/np)*np);
end

dims = size(image);
% introduce planar phase
[xx,yy] = ndgrid(-dims(1)/2:dims(1)/2-1,-dims(2)/2:dims(2)/2-1);
ph_max = pi/5;
ph = ph_max*xx/(dims(1)/2)+ph_max*yy/(dims(2)/2);
image = image.*exp(1i*ph);

% generate sensitivity maps
smap = mri_sensemap_sim('nx',dims(1),'ny',dims(2),'ncoil',nc);
close;

im_rep = repmat(image,[1 1 nc]);

% apply sensitivity maps
mapped_im = im_rep.*smap;

% throw away k-space data
im_fft = zeros(size(mapped_im));
if (reduced_dim == 1)
    reduced_fft = zeros(dims(1)/np,dims(2),nc);
else
    reduced_fft = zeros(dims(1),dims(2)/np,nc);
end
for ii = 1:nc
    im_fft(:,:,ii) = fft2(mapped_im(:,:,ii));
    if (reduced_dim == 1)
        reduced_fft(:,:,ii) = im_fft(1:np:end,:,ii);
    else
        reduced_fft(:,:,ii) = im_fft(:,1:np:end,ii);
    end
end

%% construct mask
mask = any(abs(mapped_im)>0.005,3);

%% perform SENSE reconstruction for varying levels of noise
sigmas = 0:25:350;
recon_im3 = zeros(dims(1),dims(2),length(sigmas));
error3 = zeros(length(sigmas),1);
if (comparison)
    recon_im1 = zeros(dims(1),dims(2),length(sigmas));
    error1 = zeros(length(sigmas),1);
    recon_im2 = zeros(dims(1),dims(2),length(sigmas));
    error2 = zeros(length(sigmas),1);
    recon_im4 = zeros(dims(1),dims(2),length(sigmas));
    error4 = zeros(length(sigmas),1);
end
for ii = 1:length(sigmas)
    % introduce complex AWGN
    sigma = sigmas(ii);
    noisy_reduced_fft = reduced_fft + sigma*randn(size(reduced_fft)) + 1i*sigma*randn(size(reduced_fft));
    
    % SENSE reconstruction
    recon_im3(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, reduced_dim, np, regularizer1, mask, figs);
    error3(ii) = sqrt((sum(sum((real(image)-real(recon_im3(:,:,ii))).^2))+sum(sum((imag(image)-imag(recon_im3(:,:,ii))).^2)))/prod(dims));
    if (comparison)
        recon_im1(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, reduced_dim, np, regularizer1, true(size(image)), figs);
        recon_im2(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, reduced_dim, np, regularizer2, true(size(image)), figs);
        recon_im4(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, reduced_dim, np, regularizer2, mask, figs);
        
        error1(ii) = sqrt((sum(sum((real(image)-real(recon_im1(:,:,ii))).^2))+sum(sum((imag(image)-imag(recon_im1(:,:,ii))).^2)))/prod(dims));
        error2(ii) = sqrt((sum(sum((real(image)-real(recon_im2(:,:,ii))).^2))+sum(sum((imag(image)-imag(recon_im2(:,:,ii))).^2)))/prod(dims));
        error4(ii) = sqrt((sum(sum((real(image)-real(recon_im4(:,:,ii))).^2))+sum(sum((imag(image)-imag(recon_im3(:,:,ii))).^2)))/prod(dims));
    end
    
end
%% plot error as a function of noise variance
figure; plot(sigmas.^2, error3, 'b');
if (comparison)
    hold on; plot(sigmas.^2, error1,'k');
    plot(sigmas.^2, error2, 'r');
    plot(sigmas.^2, error4, 'g');
    legend('no reg, mask', 'no reg, no mask', 'tikhonov, no mask', 'tikhonov, mask');
end
xlabel('variance of complex AWGN in kspace');
ylabel('RMS error in reconstructed images');
if (comparison)
    title('effect of Tikhonov regularization on noisy SENSE reconstruction')
else 
    title('RMSE as function of noise variance');
end
%% plot reconstructed image at sigma(ii)
ii = 6; % change as desired
if (comparison)
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
else
    figure; subplot(231); imshow(real(image),[]); colorbar;
    title('real of original image');
    subplot(232); imshow(imag(image),[]); colorbar;
    title(['imag of original image, noise variance: ' num2str(sigmas(ii).^2)]);
    subplot(234); imshow(real(recon_im3(:,:,ii)),[]); colorbar;
    title('real of recon image, masked');
    subplot(235); imshow(imag(recon_im3(:,:,ii)),[]); colorbar;
    title('imag of recon image, masked');
    subplot(236); imshow(abs(real(image)-real(recon_im3(:,:,ii)))+abs(imag(image)-imag(recon_im3(:,:,ii))),[]); colorbar;
    title('abs(differences), masked');
end
