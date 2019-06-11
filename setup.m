% setup.m
% run this file to set up matlab path etc.
% you may need to modify this depending on how you installed the toolbox
% so this should be considered simply a "guide" not a robust script.

if ~exist('irtdir', 'var')
	disp('The variable "irtdir" is not set, so trying default, assuming')
	disp('that you launched matlab from the irt install directory.')
	disp('You may need to edit setup.m or adjust your path otherwise.')

%	irtdir = pwd; % the default is to assume launch from irt directory

	% default is to look for directory where this setup.m is installed!
	irtdir = which('setup'); % find setup.m
	[irtdir, ~] = fileparts(irtdir);

	disp(['Assuming you installed irt in directory "' irtdir '".'])

%	irtdir = '~fessler/l/src/matlab/alg/'; % where you install this package
%	irtdir = '~fessler/l/web/irt/'; % where you install this package
end

if ~exist(irtdir, 'dir')
	disp(sprintf('The directory "%s" does not exist', irtdir))
	error(sprintf('you need to edit %s to change default path', mfilename))
end

if irtdir(end) ~= filesep % make sure there is a '/' at end of directory
	irtdir = [irtdir filesep];
end

list = {...
'align', ...		% image registration
'align/mex', ...	% image registration mex files
'blob', ...		% blob (KB) basis
'ct', ...		% x-ray CT (polyenergetic) recon
'data', ...		% example data
'emission', ...		% emission image reconstruction
'example', ...		% example applications
'fbp', ...		% FBP (filtered backprojection) code
'general', ...		% generic image reconstruction
'graph', ...		% graphics routines
'mri', ...		% MRI reconstruction
'mri-rf/yip-spsp', ...	% MRI RF pulse design
'nufft', ...		% nonuniform FFT (for a fast projector)
'nufft/table', ...	% mex files for NUFFT
'penalty', ...		% regularizing penalty functions
'systems', ...		% system "matrices"
'systems/tests', ...	% tests of systems
'transmission', ...	% transmission image reconstruction
'utilities', ...	% various utility functions
'wls', ...		% weighted least-squares (WLS) estimates
'um', ...		% extra directories for UM users only
'um/lustig-poisson-disk'
};
%'mri/recon', ...	% MRI reconstruction - old

for ii=1:numel(list)
	tmp = [irtdir list{ii}];
	if exist(tmp, 'dir'), addpath(tmp), end
end


% Set up path to mex files, possibly depending on matlab version.
% Fortunately it seems that the v6 files also run on v7.
% Unfortunately, it seems that sometimes compiling on one version
% e.g., 7.3 and running on earlier 7.0 won't work.
% If you have mex problems, comment out the following line.
% Much of the toolbox will work without mex, just slower.
addpath([irtdir 'mex/v7']);


% do not add the paths below if you are using Matlab!
if ir_is_octave
	addpath([irtdir 'octave']); % extra stuff for octave only!
elseif isempty(which('dbstack')) % for freemat only!
	addpath([irtdir 'freemat']); % extra stuff for freemat only!
end

% check to see if path setup worked by looking for im() routine.
if strcmp([irtdir 'graph' filesep 'im.m'], which('im'))
	disp('Path setup for irt appears to have succeeded.')
	clear list ii irtdir tmp
else
	disp('Path setup for irt may have failed.')
end
