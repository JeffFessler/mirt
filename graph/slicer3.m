 function slicer3(vol, varargin)
%function slicer3(vol, varargin)
% show transaxial, coronal, and saggital slices of a 3D volume
% in
%	vol	[nx,ny,nz]		3d image volume
% option
%	'how'	'c|cent|center|max'
%	'args'	{}			for im(..., args{:})
%	'clim'	[2]			color limits
%
% Copyright 2007-1-28, Jeff Fessler, The University of Michigan
if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(vol, 'test'), slicer3_test, return, end

arg.how = 'center';
arg.clim = [];
arg.args = {'cbar'};
arg = vararg_pair(arg, varargin);

[nx ny nz] = size(vol);

switch arg.how
case {'c', 'cent', 'center'}
	% position of center point (where 3 planes intersect)
	cx = round((nx+1)/2);
	cy = round((ny+1)/2);
	cz = round((nz+1)/2);
case 'max'
	tmp = imax(vol, 3);
	[cx cy cz] = deal(tmp(1), tmp(2), tmp(3));
otherwise
	fail('bad arg.how "%s"', arg.how)
end

if ~isempty(arg.clim)
	arg.args = {arg.args{:}, arg.clim};
end

im clf
im pl 2 2
im(1, vol(:,:,cz), sprintf('Transaxial, cz=%g', cz), arg.args{:})
im(2, vol(:,cy,:), sprintf('Coronal?, cy=%g', cy), arg.args{:})
im(3, vol(cx,:,:), sprintf('Sagital?, cx=%g', cx), arg.args{:})


%
% slicer3_test
%
function slicer3_test
ig = image_geom('nx', 2^6, 'ny', 2^6-2, 'nz', 2^5, 'fov', 240, 'zfov', 90);
vol = ellipsoid_im(ig, [], 'oversample', 2);
clim = [0.98 1.05];
% clf, im(vol, clim), cbar
slicer3(vol, 'args', {'cbar'}, 'clim', clim, 'how', 'c')
im(4, vol, clim), cbar
