function recon_im = SENSE_recon(smap, reduced_fft, reduced_dim, np, regularizer, figs_on)
% function recon_im = SENSE_recon(smap, reduced_fft, reduced_dim, np,
%       regularizer, figs_on)
% in:
%   smap          [N M nc]                complex sensitivity maps
%   reduced_fft   [N/np M] or [N M/np]    undersampled fft
%   reduced_dim   1 or 2                  dimension of undersampling
%   np            degree of undersampling
%   regularizer   string, 'none' or 'tikhonov'
%   figs_on       boolean that shows intermediate images for debugging
% out:
%   recon_im      [N M]                   SENSE reconstructed image
% This implementation assumes that the fft will be reduced in such a way
% that the DC term is still included (i.e. even DFT coefficients retained)
% 2012-06-07 Mai Le, University of Michigan

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
beta = 0.01*max(SIG(:));
beta = 0.005*max(SIG(:));

% reconstruct image by finding (regularized) LS solution for each set of pixels
recon_im = zeros(dims);
if (reduced_dim == 1)
    for ii = 1:dims(1)/np
        for jj = 1:dims(2)
            S = constructS(reduced_dim, smap, ii, jj, np);
            a = squeeze([aliased_im(ii,jj,1:nc)]);
            v = recon_pixels(S,a,regularizer, beta);
            recon_im = place_pixels(recon_im, v, reduced_dim, dims, ii, jj, np);
        end
    end
else
    for ii = 1:dims(1)
        for jj = 1:dims(2)/np
            S = constructS(reduced_dim, smap, ii, jj, np);
            a = squeeze([aliased_im(ii,jj,1:nc)]);
            v = recon_pixels(S,a,regularizer, beta);
            recon_im = place_pixels(recon_im, v, reduced_dim, dims, ii, jj, np);
        end
    end
end

% selects values from sensitivity map relevant for specific pixels
function S = constructS(reduced_dim, smap, ii, jj, np)
dims = size(smap);
indeces = (0:np-1)*dims(reduced_dim)/np;
if (reduced_dim == 1)
    S = squeeze([smap(ii+indeces,jj,:)]).';
else
    S = squeeze([smap(ii,jj+indeces,:)]).';
end

% places each of the reconstructed pixels in the reconstructed image
function recon_im = place_pixels(recon_im, v, reduced_dim, dims, ii, jj, np)
indeces = (0:np-1)*dims(reduced_dim)/np;
if (reduced_dim == 1)
    for kk = 1:np
        recon_im(ii+indeces(kk),jj) = v(kk);
    end
else
    for kk = 1:np
        recon_im(ii,jj+indeces(kk)) = v(kk);
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

