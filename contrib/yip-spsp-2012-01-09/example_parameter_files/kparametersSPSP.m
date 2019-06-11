 function [kp] = kparametersSPSP()
%function [kp] = kparametersSPSP() kp is a structure that contains
%excitation k-space trajectory parameters for spectral-spatial pulse
%design. The trajectory comes from an oscillatory Gz waveform comprised of
%trapazoids with opposite polarities, which is most commonly used for 
%conventional SPSP pulses.
%
% Note: default parameter values as on our paper. They will generate
% the alpha=-2E-4 pulse, as shown in Fig. 3b.
%
%Chun-yu Yip, 4/1/2009

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