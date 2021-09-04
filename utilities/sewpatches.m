function [ output_im ] = sewpatches(input_patches, output_sz, patch_stsz)

patch_sz = size(input_patches,1);
[~, patch_indices] = createpatches(zeros(output_sz), patch_sz, patch_stsz);

npatches_per_image = size(patch_indices,3);
npatches = size(input_patches,3);
nrow = output_sz(1);
ncol = output_sz(2);

nims = npatches / npatches_per_image;

output_im = zeros(nrow,ncol,nims);

for ii = 1:nims
    psf_img = zeros(nrow, ncol);
    sew_img = zeros(nrow, ncol);
    for indx = 1:npatches_per_image
        patch = reshape(input_patches(:,:,(ii-1)*npatches_per_image + indx),[patch_sz, patch_sz]);
        
        ind = patch_indices(:,:,indx);
        r1 = ind(1); c1 = ind(2); r2 = ind(3); c2 = ind(4);
        
        psf_img(r1:r2,c1:c2) = psf_img(r1:r2,c1:c2) + 1;
        sew_img(r1:r2,c1:c2) = sew_img(r1:r2,c1:c2) + patch;
        
    end
    output_im(:,:,ii) = sew_img ./ (psf_img + 1e-12);
end


end