function [iop] = ioparametersSPSP()
%function [iop] = ioparametersSPSP()
%rfp is a structure that contains parameters for input and 
%output of the pulse design code.
%
%Chun-yu Yip, 4/1/2009

iop.writetofile_sim = logical(1);       %Write designed pulse and gz waveforms
                                        %to files for Bloch simulation.
                                                                             
iop.writetofile_scanner = logical(1);   %Write designed pulse and gz waveforms
                                        %to files for MRI scanner.

%Bloch simulation parameters
iop.bfovz = 3;                          %cm; range of z (space) in Bloch
                                        %simulation = [-bfovz/2, bfovz/2).
iop.bdimz = 200;                        %Number of z samples in Bloch simulation.
iop.bfovf = 500;                        %Hz; range of f (frequency) in Bloch
                                        %simulation = [-bfovf/2, bfovf/2).
iop.bdimf = 200;                        %Number of f samples in Bloch simulation.

%File names
iop.infname = 'spspexample';            %File name for RF and z gradient waveform 
                                        %files. Extension should not be included.
iop.outfname_from_sim = 'mspsp';        %File name for Bloch simulation results.
                                        %Extension should not be included.

iop.show_blochsimresults = logical(1);  %Display Bloch simulation results

%Where pulse gradient waveforms and simulation results reside. 
cd('..');
[status,homepath] = system('pwd');
homepath = strtrim(homepath);
iop.waveformfilespath = [homepath '/example_waveform_files/'];

