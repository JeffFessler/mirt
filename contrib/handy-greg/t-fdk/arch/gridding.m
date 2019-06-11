function [newProj, rectGrid] = gridding(cg, proj )

indexT = 1:cg.nt;
indexS = 1:cg.ns;

ws = (cg.ns+1)/2;
wt = (cg.nt+1)/2;

coordT = (indexT - wt) * cg.dt;
coordS = (indexS - ws) * cg.ds;
d = cg.dsd;
r = cg.dso;

matrix = zeros(cg.nt, cg.ns);
%create matrix
for i = 1:1:cg.ns
	matrix(:,i) = coordT .* sqrt(r^2-coordS(i)^2) / d;
end

%find the maximum value
upperLimit = max(max(matrix));
lowerLimit = min(min(matrix)); % should just be -max

%create new grid in the vertical direction
new_dt = (upperLimit-lowerLimit)/cg.nt; % original
%new_dt = cg.dt;
rectGrid = lowerLimit+(upperLimit-lowerLimit)/cg.nt:new_dt:upperLimit;
rectGrid = rectGrid'; % 1d set of "t" samples we want

%must also loop over the projection angles
newProj = zeros(cg.ns, numel(rectGrid), cg.na);
for j = 1:1:cg.na
	for i = 1:1:cg.ns
		newProj(i,:,j) = ...
		interp1(matrix(:,i), proj(i,:,j), rectGrid, 'linear', 'extrap');
	end
end

end % gridding()
