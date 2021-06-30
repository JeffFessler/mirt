% ir_mex_build
%
% Some of the functions in this toolbox use calls to mex files.
% The source code for *some*, but not all, of those mex files
% is included with the toolbox, and this script generates the corresponding mex files.

% The full MIRT toolbox with complete mex files can be found at:
% http://web.eecs.umich.edu/~fessler/code/index.html

% This script and the scripts it calls have been tested
% for compiling on Mac, PC, and Linux.

ir_mex_build_mri
ir_mex_build_table

% todo: ir_mex_build_penal
