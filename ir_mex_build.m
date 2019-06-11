% ir_mex_build
%
% some of the functions in this toolbox use calls to mex files.
% the source code for *some*, but not all, of those mex files
% is included with the toolbox.
% the compiled mex files are available for several supported system types
% in the mex directory.  for unsupported systems, e.g., PCs, the user
% will have to compile her own mex files using this script.
%
% This script and the scripts it calls are largely untested
% because I use the Makefiles in each directory
% instead of the Matlab command-line "mex" command
% for compiling on Mac and Linux.

ir_mex_build_mri
ir_mex_build_table

% todo: ir_mex_build_penal
