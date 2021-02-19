% ir_mex_build_mri
% run matlab's "mex" command to "compile" the mri-related code
% into mex files.
% only users on unsupported systems, e.g., PCs, will need to do this

% I never use this file because I use "make" with ./Makefile instead.
% This is just an example.
% You might have to modify the flags for your system.

dir_current = pwd;
dir_mri = path_find_dir('mri');
cd(dir_mri)

% compile mex files using ir_mex_fun
fun = ir_mex_fun();
fun('exp_xform_mex.c')
fun('mri_exp_mult_mex.c')

cd(dir_current)

% "test" mex files by running them (each should display usage)
try 
	exp_xform_mex
	mri_exp_mult_mex
catch
    disp(['WARNING! Error compiling mex files in ', mfilename])
end
