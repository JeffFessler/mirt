function [patch_ims, patch_indices] = createpatches(full_im, patch_sz, patch_stsz)


[nrow, ncol, nims] = size(full_im);
patch_hs = (patch_sz - 1) / 2;  % patch half size
rindx = patch_hs:patch_stsz:nrow-patch_hs-1;
cindx = patch_hs:patch_stsz:ncol-patch_hs-1;
npatch = nims * numel(rindx) * numel(cindx);
patch_ims = zeros(patch_sz, patch_sz, npatch);
patch_indices = zeros(2,2,npatch);

indx = 1;
for nn = 1:nims
    for x = 1:numel(rindx)
        for y = 1:numel(cindx)
            r1 = rindx(x)+1-patch_hs; r2 =  rindx(x) + patch_hs + 1;
            c1 = cindx(y)+1-patch_hs; c2 = cindx(y) + patch_hs + 1;
            patch = full_im(r1:r2, c1:c2, nn);
            patch_ims(:, :, indx) = patch;
            patch_indices(:,:,indx) = [r1,r2;c1,c2];
            indx = indx + 1;
        end
    end
end


end