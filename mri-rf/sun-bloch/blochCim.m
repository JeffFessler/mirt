%	Bloch simulator using mex file
%   [mx,my,mz] = blochCim(b1,gr,tp,T1,T2,freq,pr,mode,sens,mx0,my0,mz0);
%
%	Bloch simulation of rotations due to B1, gradient and
%	off-resonance, including relaxation effects.  At each time
%	point, the rotation matrix and decay matrix are calculated.
%	Simulation can simulate the steady-state if the sequence
%	is applied repeatedly.This code can also simulate the parallel transmit.
%
%   INPUT:
%   b1 = (ntime * ncoil) RF pulse, complex number in (Gauss).
%   gr = (ntime * 3) gradient in (Gauss/cm).
%   tp = (ntime * 1) time duration of each b1 and gr point, in (sec),
%				or (1 * 1) time step if constant for all points
%				or monotonically INCREASING endtime of each
%				interval..
%   T1, T2 = (1 * 1) relaxation time in (sec)
%   freq = (nx * ny * nz) or (npos * 1) off-resonance freq. in (Hz)
%   pr = (npos * 3) array of spatial positions in (cm)
%   mode= Bitmask simulation mode: (default 0)
%		Bit 0:  0-Simulate from start, 1-Steady State
%		Bit 1:  1-Record m at time points.  0-just end time.%
%       Bit 2:  0-display simulation information. 1-only print out warnings.
%       
%
%   (OPTIONAL)
%   sens = (ncoils * npos or ncoils * nx*ny*nz) sensitivities matrix.
%   default is 1 for all position. Set to default if you put 1 here.  
%   mx0,my0,mz0 = (npos * 1) initial magnetization
%   default is (0,0,1) for all position
%
%	OUTPUT:
%		mx,my,mz = npos x 1 arrays of the resulting magnetization
%				components at each position.
%
%	The code is originally downloaded from Brian Hargreaves's website.
%   
%   Hao Sun made the following main changes:
%   (1)match field map with position
%   (2)correct the minus sign error in the 'b1imag' in the blochsimfz sub
%   function in Brian's code
%   (3)add parallel simulation function
%   
%
%   Hao Sun, the University of Michigan, Jul 22,2012






