function recon_im = sense_recon(smap, reduced_fft, reduced_dim, np, regularizer, mask, figs_on)
% function recon_im = sense_recon(smap, reduced_fft, reduced_dim, np,
%       regularizer, figs_on)
% in:
%   smap          [N M nc]                complex sensitivity maps
%   reduced_fft   [N/np M] or [N M/np]    undersampled fft
%   reduced_dim   1 or 2                  dimension of undersampling
%   np            degree of undersampling
%   regularizer   string, 'none' or 'tikhonov'
%   mask          [N M] boolean
%   figs_on       boolean that shows intermediate images for debugging
% out:
%   recon_im      [N M]                   SENSE reconstructed image
% This implementation assumes that the fft will be reduced in such a way
% that the DC term is still included (i.e. even DFT coefficients retained)
% 2012-06-07 Mai Le, University of Michigan

assert(isequal(size(mask), size(squeeze(smap(:,:,1)))), 'mask size does not match image size');

% dimensions of original image
dims = [size(smap,1) size(smap,2)];

% number of coils
nc = size(smap,3);

% create aliased intermediate images
aliased_im = zeros(size(reduced_fft,1), size(reduced_fft,2), nc);
for ii = 1:nc
     aliased_im(:,:,ii) = ifft2(reduced_fft(:,:,ii));
end

% plot aliased images
if (figs_on)
    figure;
    for ii = 1:nc
        subplot(2,2,ii); imshow(abs(aliased_im(:,:,ii)),[]); colorbar;
    end
end

% calculate regularization parameter
% I found empirically that somewhere between 0.5% and 1% of the max SVD
% value works well for beta
C = Cdiff1(np);
[U, SIG, V] = svd(C'*C);
% beta = 0.01*max(SIG(:));
beta = 0.005*max(SIG(:));

% reconstruct image by finding (regularized) LS solution for each set of pixels
recon_im = zeros(dims);
indeces = (0:np-1)*dims(reduced_dim)/np;
if (reduced_dim == 1)
    for ii = 1:dims(1)/np
        for jj = 1:dims(2)
            if any(mask(ii+indeces,jj))
%                 S = constructS(reduced_dim, smap, ii+indeces, jj, np);
                ivect = ii+indeces;
                S = squeeze([smap(ivect(mask(ii+indeces,jj)),jj,:)]).';
                if ~all(mask(ii+indeces,jj))
                    S = S.';
                end
                a = squeeze([aliased_im(ii,jj,1:nc)]);
                v = recon_pixels(S,a,regularizer, beta);
                [recon_im] = place_pixels(recon_im, v, reduced_dim, dims, ii, jj, np, mask(ii+indeces,jj));
            end
        end
    end
else
    for ii = 1:dims(1)
        for jj = 1:dims(2)/np
            if any(mask(ii,jj+indeces))
%                 S = constructS(reduced_dim, smap, ii, jj+indeces, np, mask(ii,jj+indeces));
                jvect = jj+indeces;
                S = squeeze([smap(ii,jvect(mask(ii,jj+indeces)),:)]).';
                if ~all(mask(ii,jj+indeces))
                    S = S.';
                end
                a = squeeze([aliased_im(ii,jj,1:nc)]);
                v = recon_pixels(S,a,regularizer, beta);
                [recon_im] = place_pixels(recon_im, v, reduced_dim, dims, ii, jj, np, mask(ii,jj+indeces));
            end
        end
    end
end
% figure; subplot(121); imshow(test_fill,[]); colorbar;
% subplot(122); imshow(real(recon_im),[]); colorbar;
% keyboard;

% selects values from sensitivity map relevant for specific pixels
% depending on reduced_dim, either ivect or jvect is a 1x2 vector
% function S = constructS(reduced_dim, smap, ivect, jvect, np, mask_vector)
% 
% 
% S = squeeze([smap(ivect,jvect(mask_vector),:)]).';

% places each of the reconstructed pixels in the reconstructed image
function [recon_im] = place_pixels(recon_im, v, reduced_dim, dims, ii, jj, np, mask)
indeces = (0:np-1)*dims(reduced_dim)/np;
if (reduced_dim == 1)
    mk = 1;
    vk = 1;
    while (mk <= np)
        if (mask(mk))
            recon_im(ii+indeces(mk),jj) = v(vk);
            vk = vk + 1;
        end
        mk = mk + 1;
    end
else
    vk = 1;
    mk = 1;
    while (mk <= np)
        if (mask(mk))
            recon_im(ii,jj+indeces(mk)) = v(vk);
            vk = vk + 1;
        end
        mk = mk + 1;
    end
end

% reconstruct pixels with or without Tikhonov regularization
function pixels = recon_pixels(S,a,regularizer, beta)
switch regularizer
    case 'none'
        pixels = S\a;
    case 'tikhonov'
        M = min(size(S));
        pixels = inv(S'*S+(beta^2)*eye(M))*S'*a;
    otherwise
        pixels = S\a;
end
