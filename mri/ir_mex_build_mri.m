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

% mex compile command with optimization and c99 flag:
fun = @(f1) mex('-O', 'CFLAGS="\$CFLAGS -std=c99"', f1);
fun('exp_xform_mex.c')
fun('mri_exp_mult_mex.c')

%{ old way
mex -O CFLAGS="\$CFLAGS -std=c99" exp_xform_mex.c
mex -O CFLAGS="\$CFLAGS -std=c99" mri_exp_mult_mex.c
%}

cd(dir_current)

if 1 % "test" by running them (each displays usage)
	exp_xform_mex
	mri_exp_mult_mex
end
