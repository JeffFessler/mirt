% getFWHM_uhe(r)
%
%   Get FWHM (mm) of the ultra-high-energy collimator PSF
%   at some distances 'r' (mm) from the collimator
%
% Qiang (Victor) Lin
% Modified by A. Yendiki, 2/5/02

function FWHM = getFWHM_uhe(r)

r(r < 0) = 0;

%  Radius (mm) FWHM (mm?)
%   20          5.6135
%   75         13.8277 
%   130        18.2010
%   185        22.2908
%   245        26.9498

fw20  = 5.6135;
fw75  = 13.8277;
fw130 = 18.2010;
fw185 = 22.2908;
fw245 = 26.9498;


FWHM = (fw75 - (75-r) * ((fw75-fw20)/(75-20))) .* (r <= 75) ...
     + (fw130 - (130-r) * ((fw130-fw75)/(130-75))) .* (r <= 130 & r > 75) ...
     + (fw185 - (185-r) * ((fw185-fw130)/(185-130))) .* (r <= 185 & r > 130) ...
     + (fw185 - (185-r) * ((fw245-fw185)/(245-185))) .* (r > 185);

