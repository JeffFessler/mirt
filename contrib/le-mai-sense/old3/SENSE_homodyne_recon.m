function [demod_recon_im, recon_im, lp_im] = SENSE_homodyne_recon(reduced_ffts, smap, overlap, ...
    homodyne_direction, SENSE_direction, regularizer, nc, np, red_dims, mask, figs)
% function [recon_im, demod_recon_im, lp_im] = SENSE_homodyne_recon(reduced_ffts, smap, overlap, ...
%    homodyne_direction, SENSE_direction, regularizer, nc, np, figs)
% combines SENSE and homodyne imaging, can handle all four cases of
% combinations undersampling directions
% inputs:
%   reduced_ffts        [N/np/2+overlap+1 M] or [N/np M/2+overlap+1] or 
%                       [N/2+overlap+1 M/np] or [N M/np/2+overlap+1]
%                           SENSE and homodyne reduced kspace
%   smap                [N M nc]        complex sensitivity maps
%   overlap             width of additional strip of pixels kept beyond DC,
%                       used in low passed image for phase correction
%   homodyne_direction  direction of homodyne undersampling
%   SENSE_direction     direction of SENSE undersampling 
%   regularizer         string, 'none' or 'tikhonov'
%   nc                  number of coils
%   np                  degree of undersampling
%   red_dims            2x1 vector, size of fft after SENSE reduction only 
%   figs                boolean that shows intermediate images for debugging
%
% outputs:
%   demod_recon_im      [N M]   reconstructed image after demodulation by
%                               reference image
%   recon_im            [N M]   reconstructed from partial kspace data, does not use low pass
%                               reference image to demodulate, used for
%                               debugging
%   lp_im               [N M]   low pass reference image used to demodulate reconstructed image
%                               output for debugging
% 2012-06-17 Mai Le, University of MIchigan


% homodyne reconstruction for 
% - reconstructed image (from merging filter and partial kspace)
% - "low pass" reference image
aliased_lp_ffts = zeros(red_dims(1),red_dims(2),nc);
aliased_recon_ffts = zeros(red_dims(1),red_dims(2),nc);
for ii=1:nc
    [recon_ph_demod, recon_im, lp_im, full_kspace]...
        = homodyne_recon(reduced_ffts(:,:,ii), red_dims(1), red_dims(2), overlap, homodyne_direction);
    aliased_lp_ffts(:,:,ii) = fft2(lp_im);
    aliased_recon_ffts(:,:,ii) = fft2(recon_im);
end

% SENSE recon of recon image and reference image
recon_im = SENSE_recon(smap, aliased_recon_ffts, SENSE_direction, np, regularizer, mask, figs);
lp_im = SENSE_recon(smap, aliased_lp_ffts, SENSE_direction, np, regularizer, mask, figs);

% phase demodulation
demod_recon_im = real(recon_im.*exp(-1i*angle(lp_im)));

end
