 function [uu vv azim polar] = ir_coord_cb_flat_to_par(ss, tt, dso, dod)
%function [uu vv azim polar] = ir_coord_cb_flat_to_par(ss, tt, dso, dod)
%|
%| convert from cone-beam flat panel to parallel-beam coordinates.
%| ss and tt must have same size; dso and dsd are scalars
%|
%| out
%|	uu,vv		transaxial and axial parallel-beam detector coordinates
%|	azim		transaxial or azimuthal angle (radians)
%|	polar		polar angle (radians)
%|
%| 2013-05-24, Jeff Fessler, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

Ds = dso;
Dd = dod;
Dc = Ds + Dd;
% dsd = dso + dod;

% pvar = ss * Ds / Dc; % kak eq 155
% zvar = tt * Ds / Dc;
% uu = pvar * Ds ./ sqrt(Ds^2 + pvar.^2); % kak eq 156
% vv = zvar * Ds ./ sqrt(Ds^2 + zvar.^2); % kak eq 158

uu = Ds * ss ./ sqrt(Dc^2 + ss.^2);
vv = Ds * tt ./ sqrt(Dc^2 + ss.^2 + tt.^2) .* Dc ./ sqrt(Dc^2 + ss.^2);

%polar = atan(zvar/Ds); % kak eq 159
polar = -atan(tt ./ sqrt(Dc^2 + ss.^2)); % trick: empirical negative

%azim = atan(pvar/Ds); % kak eq 158
azim = atan(ss / Dc);

end % ir_coord_cb_flat_to_par()
