 function sino = sino_mash(sino, nr, nv, orbit)
%function sino = sino_mash(sino, nr, nv, orbit)
%	nr	radial sinogram mash factor
%	nv	angular view sinogram mash factor

if ~isvar('nr'), nr = 1; end
if ~isvar('nv'), nv = 1; end
if ~isvar('orbit'), orbit = 360, end

if orbit ~= 360, error 'non-360 not done', end

%
%	radial mashing - average pairs of neighboring radial elements
%
nr = round(log(nr) / log(2));
for ii=1:nr
	if rem(size(sino,1),2), error 'need even # of radial bins', end
	sino = (sino(1:2:end,:) + sino(2:2:end,:)) / 2;
end

%
%	angular mashing
%	trick: use 1/2 of the adjacent sinogram rows to preserve 0 degrees
%
nv = round(log(nv) / log(2));
for ii=1:nv
	if rem(size(sino,2),2), error 'need even # of radial bins', end
	sino = sino(:,1:2:end)/2 + sino(:,2:2:end)/4 + ...
		[sino(:,end), sino(:,2:2:end-2)]/4;
end
