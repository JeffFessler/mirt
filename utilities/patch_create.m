function [patch_ims, patch_indices] = patch_create(full_im, patch_sz, patch_stsz)
% function [patch_ims, patch_indices] = patch_create(full_im, patch_sz, patch_stsz)
%   Create smaller image patches from input images, can have a 2D image or 
%   a stack of 2D images as input.
% 
%   Note: patch_create and patch_sew differ from im2col and col2im by 
%   allowing a custom stride, and were designed to output a stack of 2d 
%   patches for direct input to a 2d neural network.
%
% Inputs:
%   full_im     [nlin, ncol, nims]  stack of input images
%   patch_sz    1x1   dimension of 2D patch (assumes square patch)
%   patch_stsz  1x1   step size btwn patches (stride)
%
% Outputs:
%   patch_ims     [nlin_patch, ncol_patch, npatch]  output patches
%   patch_indices [2x2xnpatch_per_image]            indices of output patches
%
%  Melissa Haskell, University of Michigan, 2021-09-08

%% Test routine
%  this tests both this function (patch_create), which creates patches, and 
%  also the function patch_sew, which puts them back together
if streq(full_im, 'test'), test_patch_funcs, return, end

if nargin < 3, ir_usage, end

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



function test_patch_funcs
 %% test_patch_funcs.m
%  function to test patch creation and recombination (or sewing) with
%  various image, patch, and step sizes
%
%  Melissa Haskell, University of Michigan, 2021-09-08

%% set sizes to test (good to have some even and odd of each)
n_all = 85:88;
stride = 1:5;
ps = 9:15;

nn = numel(n_all); ns = numel(stride); nps = numel(ps);
ntestcases = nn*ns*nps; ind = 0;

%% loop through various image, patch, and step sizes to test
fprintf('Testing patch creation and recombination '); nbytes = 0;
for ii = 1:nn
    for jj = 1:ns
        for kk = 1:nps
            
            ind = ind+1;
            fprintf(repmat('\b',1,nbytes))
            nbytes = fprintf('%d of %d', ind, ntestcases);

            n = n_all(ii);
            s = stride(jj);
            patch_sz = ps(kk);
            
            % create 3d test image and pad
            test_im = cat(3,phantom(n,n),phantom(n,n)*3,phantom(n,n)*8);
            test_im = padarray(test_im,[patch_sz,patch_sz],0,'both');
            n2 = size(test_im,1);
            
            % create patches and sew back together
            test_im_patches = patch_create(test_im,patch_sz,s);
            sew_test = patch_sew(test_im_patches,[n2,n2],s);
            
            % see if there is any significant error
            tmperr = norm(test_im(:)-sew_test(:))/norm(test_im(:));
            if tmperr > 1e-10
                warning(['patch test failed for image size ', num2str(n2),...
                    ', patch size ', num2str(patch_sz), ...
                    ', and stride ', num2str(s)])
            end
        end
    end
end
disp(' ')
disp('Patch test complete.')




end
