 function [kp, iop, rfp] = example_spsp_param1()
%function [kp, iop, rfp] = example_spsp_param1()
% setup parameters for spectral-spatial pulse design.
%
%Chun-yu Yip, 4/7/09

kp = kparametersSPSP();                 % k-space trajectory parameters
iop = ioparametersSPSP();               % input/ouput parameters
rfp = rfparametersSPSP();               % RF pulse design parameters


function [iop] = ioparametersSPSP()
%function [iop] = ioparametersSPSP()
%iop is a structure that contains parameters for input and
%output of the pulse design code.

iop.writetofile_sim = logical(1);       %Write designed pulse and gz waveforms
                                        %to files for Bloch simulation.

iop.writetofile_scanner = logical(1);   %Write designed pulse and gz waveforms
                                        %to files for MRI scanner.

%Bloch simulation parameters
iop.bfovz = 3;                          %cm; range of z (space) in Bloch
                                        %simulation = [-bfovz/2, bfovz/2).
iop.bdimz = 200;                        %Number of z samples in Bloch simulation
iop.bfovf = 500;                        %Hz; range of f (frequency) in Bloch
                                        %simulation = [-bfovf/2, bfovf/2).
iop.bdimf = 200;                        %Number of f samples in Bloch simulation

%File names
iop.infname = 'spspexample';            %File name for RF and z gradient waveform
                                        %files. Extension should not be included.
iop.outfname_from_sim = 'mspsp';        %File name for Bloch simulation results.
                                        %Extension should not be included.

iop.show_blochsimresults = logical(1);  %Display Bloch simulation results

iop.waveformfilespath = '.';            %Where pulse gradient waveforms and
                                        %simulation results reside.




 function [rfp] = rfparametersSPSP()
%function [rfp] = rfparametersSPSP()
%rfp is a structure that contains SPSP RF waveform computation parameters.
% Note: default parameter values as on our paper. They will generate
% the alpha=-2E-4 pulse, as shown in Fig. 3b.


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




%function [kp] = kparametersSPSP() kp is a structure that contains
%excitation k-space trajectory parameters for spectral-spatial pulse
%design. The trajectory comes from an oscillatory Gz waveform comprised of
%trapazoids with opposite polarities, which is most commonly used for
%conventional SPSP pulses.
%
% Note: default parameter values as on our paper. They will generate
% the alpha=-2E-4 pulse, as shown in Fig. 3b.

function [kp] = kparametersSPSP()
kp.pointtime = 4E-6;            %s; RF pulse sampling period,
                                %assuming it is same for z gradient.
                                %For example, for Siemens Trio, it is
                                %2.5E-6 s. For GE Signa (3T) in Michigan,
                                %it was 4E-6 s.

kp.gmax = 4;                    %g/cm; maximum gradient amplitude
kp.dgdtmax = 15000;             %g/cm/s; maximum gradient slew rate

kp.T = 0.003;                   %s; z gradient oscillation period
kp.Ntraps = 7;                  %Number of trapezoids in z gradient waveform
                                %Note:
                                %These two parameters determine the shape
                                %of the Gz waveform. Suitable values of
                                %these two values depend on 1. slice
                                %profile shape and thickness, 2. "alpha", 3.
                                %TE(and TD), and perhaps more. As I have
                                %not come up with formuli for them,
                                %I always find good values for them by trial
                                %and error (ie, put in some reasonable
                                %values, do the pulse design, and perform Bloch
                                %simulation to see if those values give me
                                %a good excitation pattern). For a certain
                                %excitation accuracy, those values should
                                %be chosen so that the pulse is as short as
                                %possible. For pulses that we show in Fig. 3
                                %of our paper, we used:
                                %   alpha=-2:   T = 0.003;    Ntraps = 7;
                                %   alpha=-2.5: T = 0.003;    Ntraps = 8;
                                %   alpha=-3:   T = 0.003;    Ntraps = 9;
                                %Other good values for different alphas:
                                %   alpha=-1:   T = 0.0025;   Ntraps = 6;
                                %   alpha=-1.25:T = 0.0025;   Ntraps = 6;
                                %   alpha=-1.5: T = 0.0025;   Ntraps = 6;
                                %   alpha=-1.75:T = 0.0025;   Ntraps = 7;
                                %   alpha=-2.25:T = 0.003;    Ntraps = 7;
                                %   alpha=-2.75:T = 0.003;    Ntraps = 8;
                                %   alpha=-3.25:T = 0.003;    Ntraps = 10;
                                %   alpha=-3.5: T = 0.003;    Ntraps = 10;

kp.pw = 9999;                   %s; pulse duration. Will be filled once the
                                %pulse is designed (9999: not known yet)
kp.npnts = 9999;                %Number of pulse samples. Will be filled once
                                %the pulse is designed (9999: not known yet)
