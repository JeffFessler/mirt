 function [uu vv azim polar] = ir_coord_cb_cp1_to_par(rr, tt, dso, dod)
%function [uu vv azim polar] = ir_coord_cb_cp1_to_par(rr, tt, dso, dod)
%|
%| convert from cone-parallel (for arc / equiangular / 3rd gen CT)
%| to parallel-beam coordinates for beta=0 which is azim=0
%| ss and tt must have same size; dso and dsd are scalars
%|
%| called "cp1" because this is the "first step" in cone-parallel rebinning
%| grass:00:3cb
%|
%| in
%|	rr	[(N)]	radial distance of rays from isocenter
%|	tt	[(N)]	axial coordinates on virtual cp1 detector
%|
%| out
%|	uu,vv		transaxial and axial parallel-beam detector coordinates
%|	azim		transaxial angle or azimuthal (radians, ccw from y-axis)
%|	polar		polar angle (radians)
%|
%| 2013-05-27, Jeff Fessler, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

uu = rr;
polar = -asin(tt / dsd);
vv = dso * sin(polar);
azim = zeros(size(rr), class(rr));

end % ir_coord_cb_cp1_to_par()
