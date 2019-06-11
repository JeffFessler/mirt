files list: 

main files:
blochCim.c  : c code for bloch simulation.
blochCim.m  : help file for blochCim.
blochCim.mexa64: compiled mex file for linux.
blochCim.mexw64: compiled mex file for Windows.

optional files: 
(demo files)
demoBlochCim.m : demo usage of main files.

(interface files)
parallel_blochCim.m: reshape and apply mask in the simulation.
steadyState_blochCim.m: simulate steady state with gradient crusher effect.

(MATLAB bloch simulation files)
compareWithMatlab.m : compare c simulation with old MATLAB code.
parallel_blochsim_field_roi.m : interface code of MATLAB bloch simulation.
blochsim2.m: core part of MATLAB bloch simulation.

You may need Jeffrey Fessler's image reconstruction toolbox to run some of the optional files. 
http://web.eecs.umich.edu/~fessler/code/index.html


Hao Sun, University of Michigan , Oct.9,2012