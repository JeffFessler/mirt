function [patch_ims, patch_indices] = createpatches(full_im, patch_sz, patch_stsz)
% function [patch_ims, patch_indices] = createpatches(full_im, patch_sz, patch_stsz)
%   Create smaller image patches from input images, can be 2D image or a
%   stack of 2D images.
%
% Inputs:
%   full_im     [nlin, ncol, nims]  stack of input patches
%   patch_sz    1x1   dimentions of 2D patch
%   patch_stsz  1x1   step size btwn patches (stride)
%
% Outputs:
%   patch_ims     [nlin_patch, ncol_patch, npatch]  output patches
%   patch_indices [2x2xnpatch_per_image]            indices of output patches
%
%  Melissa Haskell, University of Michigan, 2021-09-08

%% Get data sizes
[nrow, ncol, nims] = size(full_im);
patch_hs = (patch_sz - 1) / 2;  % patch half size

%% Initialize row and column locations for the center of the patches
rindx = patch_hs:patch_stsz:nrow-patch_hs-1;
cindx = patch_hs:patch_stsz:ncol-patch_hs-1;

%% Initialize output variables
npatch = nims * numel(rindx) * numel(cindx);
patch_ims = zeros(patch_sz, patch_sz, npatch);
patch_indices = zeros(2,2,npatch);

%% Look through each image and create patches
indx = 1;
for nn = 1:nims
    for x = 1:numel(rindx)
        for y = 1:numel(cindx)
            r1 = rindx(x)+1-patch_hs; r2 =  rindx(x) + patch_hs + 1;
            c1 = cindx(y)+1-patch_hs; c2 = cindx(y) + patch_hs + 1;
            patch = full_im(r1:r2, c1:c2, nn);
            patch_ims(:,:,indx) = patch;
            patch_indices(:,:,indx) = [r1,r2;c1,c2];
            indx = indx + 1;
        end
    end
end


end