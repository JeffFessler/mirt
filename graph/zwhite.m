 function zwhite
%function zwhite
%	in black background mode, the z axis does not show up properly
%	and this fixes it for current plot
set(gca, 'zcolor', [1 1 1])
