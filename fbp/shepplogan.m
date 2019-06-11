 function [image, ellpar] = shepplogan(nx, ny, isemis)
%function [image, ellpar] = shepplogan(nx, ny, isemis)
% convenience interface to ellipse_im.  obsolete.

if nargin < 1, nx = 64; end
if nargin < 2, ny = nx; end
if nargin < 3, isemis = 0; end

ig = image_geom('nx', nx, 'ny', ny, 'dx', 1);
if isemis
	arg = 'shepplogan-emis';
else
	arg = [];
end
[image ellpar] = ellipse_im(ig, arg);
