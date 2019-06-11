%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spectral-spatial RF pulse design script in Matlab, based on
%"Spectral-spatial RF pulse design for through-plane phase precompensatory
%slice selection for T2*-weighted functional MRI", Chun-yu Yip et al,
%Magnetic Resonance in Medicine, 2009. 
%
%Written by Chun-yu Yip, University of Michigan, Ann Arbor, 4/1/2009
%Current affiliation: Athinoula A. Martinos Center for Biomedical Imaging,
%Massachusetts General Hospital, Harvard Medical School, Charlestown, MA,
%USA. 
%Email address: chunyuy@nmr.mgh.harvard.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Running this script requires that the image reconstruction toolbox be
%installed first: http://www.eecs.umich.edu/~fessler/irt/fessler.tgz. In
%your Matlab workspace, please "addpath" the folders in the Fessler
%toolbox. For example:
%
%addpath('fesslertoolboxlocation/utilities/');
%addpath('fesslertoolboxlocation/systems/');
%addpath('fesslertoolboxlocation/nufft/');
%addpath('fesslertoolboxlocation/wls/');
%addpath('fesslertoolboxlocation/mex/');
%addpath('fesslertoolboxlocation/mex/v7/');
%etc.
%
%You can add those lines in your startup.m file so that you do not have to
%add the paths manually everytime. You can type "path" at the matlab prompt
%to check that the paths are successfully added.
%
%
%To design pulse, you have to first load pulse design parameter values. To
%start, you can adopt those we used in our paper in the
%example_parameter_files folder. To keep track of parameters, we recommend
%that you create one folder for each pulse design occasion. For example,
%you can create a "siemens_parameter_files" folder and a
%"ge_parameter_files", or "parameter_files_1Jan2010".
%

[status,homepath] = system('pwd');
homepath = strtrim(homepath);
cd([homepath '/example_parameter_files/']);

kp = kparametersSPSP;           %k-space trajectory parameters
rfp = rfparametersSPSP;         %Pulse design parameters
iop = ioparametersSPSP;         %Input-output parameters

cd(homepath);
%
%kp, rfp, iop are matlab structures, whose fields can be accessed by, e.g.,
%"kp.pointtime". Please see parameter files for details of each field.

%Design z-gradient waveform, based on parameters in kp.
[kp,gz,kz,kf] = compute_gz_spsp(kp);

%Design complex-valued RF waveform iteratively using conjugate gradient
[b] = compute_rf_spsp_mgh(kp,rfp,gz,kz,kf);

%Write computed waveforms to files for simulation and/or scanner.
write2files_spsp(kp,rfp,iop,gz,b);

%Perform Bloch simulation in SPSP space
[mresult] = dosim7_spsp(kp,rfp,iop);

