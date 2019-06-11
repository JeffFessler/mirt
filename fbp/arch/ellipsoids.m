 function [phantom, params] = ellipsoids(nx, ny, nz, params, ...
	dx, dy, dz, varargin)
%function [phantom, params] = ellipsoids(nx, ny, nz, params, ...
%	dx, dy, dz, varargin)
%
% generate ellipsoids phantom image from parameters:
% [x_center y_center z_center x_radius y_radius z_radius
%	xy_angle_degrees z_angle_degrees amplitude [oversample]]
%
% in
%	(dx,dy,dz)	voxel size
% option
%	oversample	over-sampling factor
% op ellipsoid in aspire with nint=3 is oversample=4 = 2^(3-1) here
%
% Copyright 2004-8-13, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
% and Jeff Fessler, The University of Michigan
if nargin < 1, help(mfilename), error(mfilename), end

if nargout > 1
	warning 'units of output params not finished'
end

if nargin == 1 && streq(nx, 'test')
	phantom = ellipsoids(2^5); % arbitrary choice
	im(phantom, 'Shepp Logan', [0.9 1.1]), cbar
	if ~nargout, clear phantom, end
return
end

warn 'ellipsoids is obsolete; use ellipsoid_im instead'

arg.oversample = 1;
arg = vararg_pair(arg, varargin);

if ~isvar('ny') || isempty(ny), ny = nx; end
if ~isvar('nz') || isempty(nz), nz = nx; end
if ~isvar('dx') || isempty(dx), dx = 1.; end
if ~isvar('dy') || isempty(dy), dy = dx; end
if ~isvar('dz') || isempty(dz), dz = dx; end

if ~isvar('params'), params = []; end

ig = image_geom('nx', nx, 'ny', ny, 'nz', nz, 'dx', dx, 'dy', -dy, 'dz', dz);
[phantom, params] = ellipsoid_im(ig, params, 'oversample', arg.oversample);
