 function st = image_geom_mri(varargin)
%function st = image_geom_mri(varargin)
%|
%| Same as image_geom() but default offsets are 0.5 so that
%| image pixel indices go from -N/2 to N/2-1.
%|
%| Copyright 2007-2-22, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), image_geom_mri_test, return, end

st = image_geom(...
	'offset_x', 0.5, ...
	'offset_y', 0.5, ...
	'offset_z', 0.5, ...
	varargin{:});


function image_geom_mri_test
image_geom_mri('nx', 64, 'dx', 2);
