function recon_im = SENSE_recon(smap, reduced_fft, np, varargin)
% function recon_im = SENSE_recon(smap, reduced_fft, np, varargin)
% inputs:
%   smap          [N M nc]                complex sensitivity maps
%   reduced_fft   [N/np M nc] or [N M/np nc]    undersampled fft
%   np       degree of undersampling
% optional inputs:
%   direction     1 or 2                  dimension of undersampling    
%   reg      type of regularizer, string, 'none' or 'tikhonov'
%   mask     [N M] boolean
%   figs     boolean that shows intermediate images for debugging

% out:
%   recon_im      [N M]                   SENSE reconstructed image
%
% This implementation assumes that the fft will be reduced in such a way
% that the DC term is still included (i.e. even DFT coefficients retained)
% 2012-06-07 Mai Le, University of Michigan
% 2013-01-02 Mai Le, tweaks

% default vals
arg.direction = 1;
arg.reg = 'none';
arg.mask = true(size(smap,1),size(smap,2));
arg.figs = 0;
arg = vararg_pair(arg, varargin); % replace default vals

assert(isequal(size(arg.mask), size(squeeze(smap(:,:,1)))), ...
    'mask size does not match image size');

if arg.direction == 2 % trick to switch directions
    for ii = 1:size(smap,3)
       smapT(:,:,ii) = smap(:,:,ii).'; 
       reduced_fftT(:,:,ii) = reduced_fft(:,:,ii).';
    end
    recon_im = SENSE_recon(smapT, reduced_fftT, np, 'direction', ...
        1, 'mask', arg.mask.', 'reg', arg.reg);
    recon_im = recon_im.';
    return
elseif arg.direction ~= 1
    fail('bad direction %d',arg.direction)
end

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
if (arg.figs)
    figure;
    for ii = 1:nc
        subplot(ceil(nc/2),2,ii); 
        imshow(abs(aliased_im(:,:,ii)),[]); 
        colorbar;
    end
end

% calculate regularization parameter
% I found empirically that somewhere between 0.5% and 1% of the max SVD
% value works well for beta
C = Cdiff1(np);
[U, SIG, V] = svd(C'*C);
beta = 0.01*max(SIG(:));
% beta = 0.005*max(SIG(:));

% reconstruct image by finding (regularized) LS solution for each set of pixels
recon_im = zeros(dims);
indeces = (0:np-1)*dims(arg.direction)/np;
% if (arg.direction == 1)
    for ii = 1:dims(1)/np
        for jj = 1:dims(2)
            if any(arg.mask(ii+indeces,jj))
                ivect = ii+indeces;
                S = squeeze(smap(ivect(arg.mask(ii+indeces,jj)),jj,:));
                if sum(arg.mask(ii+indeces,jj)) > 1
                    % because if S is vector, squeeze flips S
                    S = S.';
                end
                a = squeeze(aliased_im(ii,jj,1:nc));
                v = recon_pixels(S,a,arg.reg, beta);
                [recon_im] = place_pixels(recon_im, v, arg.direction, ...
                    dims, ii, jj, np, arg.mask(ii+indeces,jj));
            end
        end
    end
% else
%     for ii = 1:dims(1)
%         for jj = 1:dims(2)/np
%             if any(arg.mask(ii,jj+indeces))
%                 jvect = jj+indeces;
%                 S = squeeze(smap(ii,jvect(arg.mask(ii,jj+indeces)),:));
%                 if sum(arg.mask(ii,jj+indeces)) > 1
%                     % because if S is vector, squeeze flips S
%                     S = S.';
%                 end
%                 a = squeeze(aliased_im(ii,jj,1:nc));
%                 v = recon_pixels(S,a,arg.reg, beta);
%                 [recon_im] = place_pixels(recon_im, v, arg.direction, ...
%                     dims, ii, jj, np, arg.mask(ii,jj+indeces));
%             end
%         end
%     end
% end

% places each of the reconstructed pixels in the reconstructed image
function [recon_im] = place_pixels(recon_im, v, direction, dims, ...
    ii, jj, np, mask)
indeces = (0:np-1)*dims(direction)/np;
if (direction == 1)
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
    case 'Tikhonov'
        M = min(size(S));
        pixels = (S'*S+(beta^2)*eye(M))\(S'*a);
    otherwise
        pixels = S\a;
end

