 function sinop = fan_rebin(sinof, dl, sd, od, nr, na, dr)
%function sinop = fan_rebin(sinof, dl, sd, od, nr, na, dr)
% rebin 360 fan beam sinogram into 360 parallel beam sinogram
% in:
%	sinof [nl,nv]	fan beam sinogram
%	dl		detector fov (lateral)
%	sd		source-detector distance
%	od		object detector distance
%	nr		desired # of radial bins (parallel beam)
%	na		desired # of angles (parallel beam)
%	dr		desired parallel radial spacing (NOT DONE)
% out:
%	sinop [nr,na]	parallel beam sinogram
%
% Copyright Apr 2000, Idris Elbakri, The University of Michigan

%
% fan beam parameters
%
[nl nv] = size(sinof);

sinof = [sinof, sinof(:,1)];	% repeat 0-degree view at 360-degrees

phi_o = [0:nv]/nv*2*pi;		% 0 - 360 angles, inclusive, since 0 repeated
l_o = [-(nl-1)/2:(nl-1)/2]' * dl;
%disp(sprintf('fan beam ray spacing = %g', dfov / nl))
disp(sprintf('fan detector fov = %g', dl * nl))

alpha = atan(l_o / sd);		% within-fan angle


%
% parallel beam parameters
%
theta = [0:na-1]*2*pi/na;
rp = [-(nr-1)/2:(nr-1)/2]' * dr;

%
% desired fan beam parameters
%
alpha_w = sin(rp/(sd-od));
for ia=1:na
	phi_w(:,ia) = theta(ia) - alpha_w;
end

phi_w(phi_w < 0) = phi_w(phi_w < 0) + 2*pi;
phi_w(phi_w > 2*pi) = phi_w(phi_w > 2*pi) - 2*pi;

l_w = sd * tan(asin(rp / od));
l_w = l_w(:,ones(1,na));

sinop = interp2(l_o, phi_o, sinof', l_w, phi_w);
