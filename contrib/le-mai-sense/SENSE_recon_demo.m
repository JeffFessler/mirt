% SENSE_recon_demo.m
% - demonstrates use of SENSE_recon.m, which implements the SENSE
% reconstruction
% - figure shows how Tikhonov regularization mitigates the increase in error
% as noise variance increases
% - applies complex AWGN to kspace
% - image must be even both directions
%
% 2012-06-08, Mai Le
% 2013-01-02 Mai Le, tweaks

% comparison = 1 will demonstrate effect of masks and regularization on MSE,
% but demo takes much longer. comparison = 0 runs masked SENSE w/o
% regularization
comparison = 0;

% number of coils
nc = 4;
% factor of undersampling
np = 2;
% dimension to undersample
for reduced_dim = [1 2]
    display(['SENSE reduced in dim ' num2str(reduced_dim)]);

    % set up "true" brain image
    %f.dir = [path_find_dir('mri') '/../data/mri/'];
    %f.xtrue = [f.dir 'brainweb_t1.jpg'];
    %mag_true = double(imread(f.xtrue))'; % true image magnitude
    mag_true = 25*ellipse_im(256);
    % requires image to be multiple of np in dimension of undersampling
    % and even for generating smap
    if (reduced_dim == 1)
        mag_true = mag_true(1:floor(size(mag_true,1)/(np*2))*(np*2),:);
    elseif (reduced_dim == 2)
        mag_true = mag_true(:,1:floor(size(mag_true,2)/(np*2))*(np*2));
    else
        fail 'bug'
    end
    

    dims = size(mag_true);
    % introduce planar phase
    [xx,yy] = ndgrid(-dims(1)/2:dims(1)/2-1,-dims(2)/2:dims(2)/2-1);
    ph_max = pi/5;
    ph = ph_max*xx/(dims(1)/2)+ph_max*yy/(dims(2)/2);
    image = mag_true.*exp(1i*ph);

    % generate sensitivity maps
    smap = mri_sensemap_sim('nx',dims(1),'ny',dims(2),'ncoil',nc,'rcoil',400);
    close;

    im_rep = repmat(image,[1 1 nc]);

    % apply sensitivity maps
    mapped_im = im_rep.*smap;

    % throw away k-space data
    im_fft = zeros(size(mapped_im));
    if reduced_dim == 1
        reduced_fft = zeros(dims(1)/np,dims(2),nc);
    elseif reduced_dim == 2
        reduced_fft = zeros(dims(1),dims(2)/np,nc);
    else
        fail 'bug'
    end
    for ii = 1:nc
        im_fft(:,:,ii) = fft2(mapped_im(:,:,ii));
        if reduced_dim == 1
            reduced_fft(:,:,ii) = im_fft(1:np:end,:,ii);
        elseif reduced_dim == 2
            reduced_fft(:,:,ii) = im_fft(:,1:np:end,ii);
        else
            fail 'bug'
        end
    end

    %% construct mask
    thresh = 0.005;
    thresh_mask = any(abs(mapped_im)>thresh,3);

    %% perform SENSE reconstruction for varying levels of noise
    f.sigmas = 0:50:550; % complex AWGN sdtd dev
    clear recon
    recon(1).im = zeros(dims(1),dims(2),length(f.sigmas));
    recon(1).err = zeros(length(f.sigmas),1);
    if comparison
        recon = repmat(recon,4,1);
    end
    nrms_fun = @(x,y) 100 * nrms(x(:), y(:));
    for ii = 1:length(f.sigmas)
        % introduce complex AWGN
        sigma = f.sigmas(ii);
        noisy_reduced_fft = reduced_fft + sigma*randn(size(reduced_fft)) + ...
            1i*sigma*randn(size(reduced_fft));

        % SENSE reconstruction
        recon(1).im(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, ...
            np, 'mask', thresh_mask, 'direction', reduced_dim);
        if comparison
            recon(2).im(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, ...
                np, 'direction', reduced_dim);
            recon(3).im(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, ...
                np, 'reg', 'Tikhonov', 'direction', reduced_dim);
            recon(4).im(:,:,ii) = SENSE_recon(smap, noisy_reduced_fft, ...
                np, 'mask', thresh_mask, 'reg', 'Tikhonov', 'direction', reduced_dim);
        end
        for num = 1:(1+3*comparison) % do all 4 if comparison
            recon(num).err(ii) = nrms_fun(image, recon(num).im(:,:,ii));
        end
    end
    %% plot error as a function of noise variance
    figure; plot(f.sigmas.^2, recon(1).err, 'bo');
    if (comparison)
        hold on; plot(f.sigmas.^2, recon(2).err,'ko');
        plot(f.sigmas.^2, recon(3).err, 'rx');
        plot(f.sigmas.^2, recon(4).err, 'gx');
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
    ii = 4; % change as desired
    clim = minmax(mag_true)';
    if (comparison)
        im(1, real(image));
        im(2, imag(image));
        titlef('imag of original image, noise variance: %d', f.sigmas(ii).^2);
        im(5, real(recon(1).im(:,:,ii)), clim);
        title('real of recon image, no reg');
        im(6, imag(recon(1).im(:,:,ii)), clim);
        title('imag of recon image, no reg');
        im(7, abs(real(image)-real(recon(1).im(:,:,ii)))+ ...
            abs(imag(image)-imag(recon(1).im(:,:,ii))), [-9 9]);
        title('abs(differences), no reg');
        im(9, real(recon(2).im(:,:,ii)), clim); 
        title('real of recon image, Tikhonov');
        im(10, imag(recon(2).im(:,:,ii)), clim); 
        title('imag of recon image, Tikhonov');
        im(11, abs(real(image)-real(recon(2).im(:,:,ii)))+ ...
            abs(imag(image)-imag(recon(2).im(:,:,ii))), [-9 9]); 
        title('abs(differences), Tikhonov');
    else
        im(1, real(image), clim); 
        title('real of original image');
        im(2, imag(image), clim); 
        titlef('imag of original image, noise variance: %d', f.sigmas(ii).^2);
        im(5, real(recon(1).im(:,:,ii)), clim);
        title('real of recon image, masked');
        im(6, imag(recon(1).im(:,:,ii)), clim);
        title('imag of recon image, masked');
        im(7, abs(real(image)-real(recon(1).im(:,:,ii)))+ ...
            abs(imag(image)-imag(recon(1).im(:,:,ii))), [-9 9]); 
        title('abs(differences), masked');
    end

prompt
end