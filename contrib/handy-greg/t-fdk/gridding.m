% gridding.m
% originally by Greg Handy
% 2013-04-04 modified by Rebecca Malinas
% 2013-04-07 streamlined by J A Fessler
function [newProj, rectGrid] = gridding(cg, proj)

indexT = 1:cg.nt;
indexS = 1:cg.ns;

ws = (cg.ns+1)/2 + cg.offset_s; % RM
wt = (cg.nt+1)/2 + cg.offset_t;

coordT = (indexT - wt) * cg.dt;
coordS = (indexS - ws) * cg.ds;
d = cg.dsd;
r = cg.dso;

matrix = zeros(cg.nt, cg.ns);
% original sample locations along t direction
matrix = coordT' * (sqrt(r^2-coordS.^2) / d); % [nt ns]

% find the maximum value
upperLimit = max(matrix(:));
lowerLimit = min(matrix(:)); % should just be -max in usual case when offset_t = 0

% create new grid in the vertical direction
nt_new = cg.nt; % todo: should allow new nt!
dt_new = (upperLimit - lowerLimit) / (nt_new-1); % new spacing

wt_new = (nt_new + 1) / 2; % trick: no offset_t in new coordintes
rectGrid = ([1:nt_new] - wt_new) * dt_new; % [nt_new] new uniform "t" samples

% loop over each projection view
newProj = zeros(cg.ns, nt_new, cg.na);
for ia = 1:cg.na
	for is = 1:cg.ns
		newProj(is,:,ia) = ...
		interp1(matrix(:,is), proj(is,:,ia), rectGrid, 'linear', 'extrap');
	end
end

end % gridding()
