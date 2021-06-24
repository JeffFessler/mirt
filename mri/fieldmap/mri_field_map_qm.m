function [wmap, wconv] = mri_field_map_qm(yik, etime, varargin)
% This function is a wrapper for fmap_est_qm.m to make it callable in a
% similar fashion as mri_field_map_reg.m

%|	yik	[(N) nset ncoil]	complex images at different echo times
%|				(these should be undistorted)
%|				(N) can be (nx,ny) or (nx,ny,nz) etc.
%|	etime	[nset]		vector of echo times [usually in seconds]
%|
%| options
%|	see fmap_est_qm.m for options
%|
%| out
%|	wmap	[(N)]	regularized estimate of field map [rad/sec]
%|	wconv	""	conventional estimate from first two scans
%|
%| Written 2021-06-03, Melissa Haskell, University of Michigan

arg = [];
arg = vararg_pair(arg, varargin, 'allow_new', 1);

%% call function
% varargin2 = arg2varargin(arg);
func_name = 'fmap_est_qm';
[wmap, wconv] = call_fmap(yik, etime, func_name, varargin{:});



end

