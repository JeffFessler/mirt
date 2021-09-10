function [ output_im ] = patch_sew(input_patches, output_sz, patch_stsz)
% function [ output_im ] = patch_sew(input_patches, output_sz, patch_stsz)
%   Combine evenly spaced patches into a single image, averaging at voxels
%   in more than one patch. Allows for the input of patches from multiple
%   output images (assuming each image has been divided into similarly 
%   sized and spaced patches.
% 
%   Note: patch_create and patch_sew differ from im2col and col2im by 
%   allowing a custom stride, and were designed to output a stack of 2d 
%   patches for direct input to a 2d neural network.
%
% Inputs:
%   input_patches   [nlin_patch, ncol_patch, npatch]  stack of input patches
%   output size     2x1   dimensions of output 2D image (or images)
%   patch_stsz      1x1   step size btwn patches (stride)
%
% Outputs:
%   output_im       [nlin, ncol, nims]   sewed together images
%
%  Melissa Haskell, University of Michigan, 2021-09-08

if nargin < 3, ir_usage, end

%% Find patch size and determine the indices of each patch using createpatches

patch_sz = size(input_patches,1);

% input 2D image to get number of patches for 2d image of that size
[~, patch_indices] = patch_create(zeros(output_sz), patch_sz, patch_stsz);

%% Calculate total number of images and initialize output

npatch_per_image = size(patch_indices,3);
npatch = size(input_patches,3);
nrow = output_sz(1);
ncol = output_sz(2);
nims = npatch / npatch_per_image;
output_im = zeros(nrow,ncol,nims);

%% Sew all the patches together
for ii = 1:nims
    psf_img = zeros(nrow, ncol);  % for tracking how many patches are at each voxel
    sum_img = zeros(nrow, ncol);  % for adding all the patch values together
    for jj = 1:npatch_per_image
        
        patch_indx = (ii-1)*npatch_per_image + jj;
        patch = reshape(input_patches(:,:,patch_indx),[patch_sz, patch_sz]);
        
        ind_2dim = patch_indices(:,:,jj);
        r1 = ind_2dim(1); c1 = ind_2dim(2); 
        r2 = ind_2dim(3); c2 = ind_2dim(4);
        
        psf_img(r1:r2,c1:c2) = psf_img(r1:r2,c1:c2) + 1;
        sum_img(r1:r2,c1:c2) = sum_img(r1:r2,c1:c2) + patch;
        
    end
    % scale summation image by the psf image
    output_im(:,:,ii) = sum_img ./ (psf_img + 1e-12);
end


end
