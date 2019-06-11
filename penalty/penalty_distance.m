  function distance = penalty_distance(offsets, sizes)
%|function distance = penalty_distance(offsets, sizes)
%|
%| convert scalar offsets to Euclidean distances to neighbors,
%| with the exception that '0' is mapped to 1 for identity case
%|
%| in
%|	offsets		[MM]
%|	sizes		[1 ndim]
%| out
%|	distance	[MM 1]
%|
%| Copyright 2006-12-6, Jeff Fessler, University of Michigan

if nargin == 1 && streq(offsets, 'test'), penalty_distance_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

displace = penalty_displace(offsets(:), sizes);
distance = sqrt(sum(displace.^2, 2));
distance(offsets == 0) = 1; % trick: for identity case


%
% penalty_distance_test()
%
function penalty_distance_test

nx = 100; ny = 80; % 2d
[ix iy] = ndgrid(-2:2, -2:2);
offsets = col(ix + iy * nx);
dis = penalty_distance(offsets, [nx ny]);
tru = sqrt(ix(:).^2 + iy(:).^2);
tru(ix == 0 & iy == 0) = 1; % trick: to match identity case
jf_equal(dis, tru)

nx = 10; ny = 8; nz = 7; % 3d
[ix iy iz] = ndgrid(-2:2, -2:2, -2:2);
offsets = col(ix + iy * nx + iz * nx * ny);
dis = penalty_distance(offsets, [nx ny nz]);
tru = sqrt(ix(:).^2 + iy(:).^2 + iz(:).^2);
tru(ix == 0 & iy == 0 & iz == 0) = 1; % trick: to match identity case
jf_equal(dis, tru)
