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
% 2013-01-02 Mai Le, tweaks

% number of coils
nc = 8;
% factor of undersampling
np = 4;

% can undersample in either direction for either method
% values should be {1,2}
for SENSE_direction = [1 2];
    for homodyne_direction = [1 2];

        % amount of overlap for homodyne reconstruction
        overlap = 20;
        regularizer = 'none'; % can be 'tikhonov' or 'none'

        % set up "true" brain image
        %f.dir = [path_find_dir('mri') '/../data/mri/'];
        %f.xtrue = [f.dir 'brainweb_t1.jpg'];
        %mag_true = double(imread(f.xtrue))'; % true image magnitude
        mag_true = 25*ellipse_im(256);
        mag_true = mag_true(7:end-6,:);
        % requires image to be multiple of np in dimension of undersampling
        if (SENSE_direction == 1)
            mag_true = mag_true(1:floor(size(mag_true,1)/(np*2))*(np*2),:);
        elseif (SENSE_direction == 2)
            mag_true = mag_true(:,1:floor(size(mag_true,2)/(np*2))*(np*2));
        else 
            fail 'bug'
        end
        dims = size(mag_true);

        % figure;
        smap = mri_sensemap_sim('nx',dims(1),'ny',dims(2),'ncoil',nc,'rcoil',400);
        close;
        
        % simulate linear phase
        [xx,yy] = ndgrid(-dims(1)/2:dims(1)/2-1,-dims(2)/2:dims(2)/2-1);
        ph_max = pi/5;
%         ph = ph_max*xx/(dims(1)/2)+ph_max*yy/(dims(2)/2);
        ph = ph_max*sin(xx/30+yy/15);
        image = mag_true.*exp(1i*ph);

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
            elseif (SENSE_direction == 2)
                reduced_fft = im_fft(:,1:np:end,ii);
            else
                fail 'bug'
            end
            SENSE_only_ffts(:,:,ii) = reduced_fft;
            reduced_fft = fftshift(reduced_fft);
            red_dims = size(reduced_fft);
            if (homodyne_direction == 1)
                reduced_fft = reduced_fft(1:min(red_dims(1),ceil(red_dims(1)/2)+overlap+1),:);
            elseif (homodyne_direction == 2)
                reduced_fft = reduced_fft(:,1:min(red_dims(2),ceil(red_dims(2)/2)+overlap+1));
            else
                fail 'bug'
            end
            reduced_ffts(:,:,ii) = reduced_fft;
        end

        %% % introduce complex Gaussian noise
        sigma = 150;
        noisy_reduced_ffts = reduced_ffts + sigma*randn(size(reduced_ffts)) + ...
            1i*sigma*randn(size(reduced_ffts));
        SENSE_only_ffts = SENSE_only_ffts + sigma*randn(size(SENSE_only_ffts))+ ...
            1i*sigma*randn(size(SENSE_only_ffts));

        %% construct mask
        mask = any(abs(mapped_im)>0.005,3);
        % patch hole in phantom which causes discontinuities in phase
        mask(60:end-60,60:end-60) = true; 

        %%
        [demod_recon_im, recon_im, lp_im] = SENSE_homodyne_recon(noisy_reduced_ffts, ...
            smap, overlap, homodyne_direction, SENSE_direction, regularizer, ...
            nc, np, red_dims, mask);
        SENSE_only_recon = SENSE_recon(smap, SENSE_only_ffts, np, ...
            'direction', SENSE_direction, 'mask', mask, 'reg', regularizer);

        %% plot SENSE & homodyne reconstructed images

        clim = minmax(mag_true)';
        im(1, mag_true, clim); 
        title('original image');

        im(2, abs(SENSE_only_recon), clim); 
        title('abs(SENSE reconstructed image)');

        im(3, mag_true-abs(SENSE_only_recon), [-50 50]);
        SENSE_only_RMS = nrms(mag_true(:),SENSE_only_recon(:));
        titlef('diff b/n SENSE only and orig, NRMS: %d', SENSE_only_RMS);

        im(5, angle(recon_im), [-3 3]);
        title('phase of recon w/o demodulation');

        im(6, real(recon_im), clim); 
        title('real of recon w/o demodulation');

        im(7, mag_true-real(recon_im), [-50 50]); 
        nondemod_RMS = nrms(mag_true(:),recon_im(:));
        titlef('diff b/n nondemod SENSE&h and orig, NRMS: %d', nondemod_RMS);

        im(9, angle(lp_im), [-3 3]); 
        title('phase of reference image');

        im(10, demod_recon_im, clim); 
        title('SENSE & homodyne, demod reconstructed image');

        im(11, mag_true-demod_recon_im, [-50 50]); 
        demod_RMS = nrms(mag_true(:),demod_recon_im(:));
        titlef('diff b/n demod SENSE&h and orig, NRMS: %d', demod_RMS);
        
        snr = sqrt(mean(abs(reduced_ffts(:)).^2))/ ...
            sqrt(mean(abs(noisy_reduced_ffts(:)).^2) - mean(abs(reduced_ffts(:)).^2));
        pr snr
        
        prompt;
    end
end
