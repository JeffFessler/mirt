function [rfp] =rfparametersSPSP()
%function [rfp] = kparametersSPSP()
%rfp is a structure that contains SPSP RF waveform computation parameters.
%
% Note: default parameter values as on our paper. They will generate
% the alpha=-2E-4 pulse, as shown in Fig. 3b.
%
%Chun-yu Yip, 4/1/2009
    
    %Slice selection parameters
    rfp.slthickz = .5 ;              %cm; slice thickness. Full width 
                                     %half max (FWHM) when shape is
                                     %"gaussian".
    rfp.flipangle = 30;              %degrees; flip angle.
    rfp.profileshape = 'rect';       %'rect' (recommended) or 'gaussian'
    rfp.smoothprofile = logical(0);  %Smoothing of the slice profile applied
                                     %to the desired pattern (not recommended).
    rfp.kernelstdv = 0.3;            %cm; Gaussian smoothing kernel std dev.
                                     %This parameter is useless when
                                     %rfp.smoothprofile = logical(0).
    rfp.hamming = logical(1);        %Hamming window in k-space applied 
                                     %directly to the designed pulse, as
                                     %described in paper. It is intended to
                                     %smooth the excited slice profile.
                                     %Recommended, although theoretically
                                     %suboptimal.
    
    %Desired pattern parameters 
    rfp.dfovz = 20;                  %cm; range of desired pattern definition
                                     %along z (space). d defined over [-dfovz/2, 
                                     %dfovz/2).
    rfp.ddimz = 1200;                %Number of samples along z.
    rfp.dfovf = 500;                 %Hz; range of desired pattern definition 
                                     %along f (frequency). f defined over
                                     %[-dfovf/2, dfovf/2).
    rfp.ddimf = 300;                 %Number of samples along f.

    %Phase pattern parameters
    rfp.zshift = -rfp.slthickz/2;    %cm; z-direction shift applied to phase 
                                     %pattern. Refer to Eq. 5 in paper.
    rfp.alpha = -2E-4;               %g/cm/Hz; key proportionality constant
                                     %between field offset and field gradient
                                     %Refer to Eq. 3 in paper.
    rfp.TE = 30E-3;                  %s; echo time of GRE sequence for fMRI.
                                     %Note that it is defined as from SPSP
                                     %pulse *center* to the time when the DC
                                     %sample is acquired.

    %Conjugate gradient parameter
    rfp.Niter = 100;                 %Number of CG iterations in design process
    rfp.beta = 0;                    %Tikhonov regularization parameter; controls 
                                     %pulse energy.