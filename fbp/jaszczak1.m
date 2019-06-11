 function ell = jaszczak1(diam)
%function ell = jaszczak1(diam)
% Generate ellipse parameters for a Jaszczak phantom of diameter diam.
% Copyright 2005-8-26, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(diam, 'test'), jaszczak1_test, return, end

nrow = [8 6:-1:2]; % how many rows of disks for each wedge
% disk radii in mm, based on Stayman, for 230mm diam disk
rads = [6.4 9 10.25 12.8 17.9 25.6] / 2 * diam / 230;

ell = [0 0 diam/2 diam/2 0 1];
for iw=1:6 % loop over six wedges
	rad = rads(iw);
	xx = []; yy = [];
	for ir=1:nrow(iw)
		x = [1:ir]' - (ir+1)/2;
		y = ir * ones(ir,1);
		xx = [xx; x];
		yy = [yy; y];
	end
	yy = yy * sqrt(3)/2;
	xx = xx * 4 * rad;
	yy = yy * 4 * rad;
	ang = iw * pi/3;
	xc = xx * cos(ang) + yy * sin(ang);
	yc = yy * cos(ang) - xx * sin(ang);
	nc = length(xc);
	ell = [ell; [xc yc rad*ones(nc,2) zeros(nc,1) -ones(nc,1)]];
end

function jaszczak1_test
ell = jaszczak1(230);
ig = image_geom('nx', 256, 'dx', 1);
x = ellipse_im(ig, ell, 'oversample', 3);
if im
	plot(ell(:,1), ell(:,2), 'o'), axis equal
	im(x)
end
