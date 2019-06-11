 function [mag, ang] = mag_angle_real(x, cutoff)
%function [mag, ang] = mag_angle_real(x, cutoff)
%
% compute the "magnitude" and "phase" of a complex number
% that is "almost" real.  so if the phase is "close" to +/- pi,
% then reverse the sign of the "magnitude" and adjust phase by pi.
% this can make for better displays of "magnitude" and "phase"
% the optional "cutoff" should be between pi/2 and pi and quantifies
% how "close" to +/- pi should be sign flipped.
%
% 2003-10-17, Jeff Fessler, The University of Michigan

if nargin < 1, ir_usage, end
if ~isvar('cutoff') || isempty(cutoff)
	cutoff = 0.6*pi;	% this "close" to +/- pi
end

mag = abs(x);
ang = angle(x);

ii = ang > cutoff;
mag(ii) = -mag(ii);
ang(ii) = ang(ii) - pi;

ii = ang < -cutoff;
mag(ii) = -mag(ii);
ang(ii) = ang(ii) + pi;
