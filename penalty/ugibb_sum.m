  function total = ugibb_sum(ugibb)
%|function total = ugibb_sum(ugibb)
%| check total wjk sum for each pixel given a ugibb stack
if nargin < 1, help(mfilename), error(mfilename), end
if streq(ugibb, 'test'), ugibb_form('test'), return, end

%
% 3D case - new way
%
if ndims(ugibb) == 4 && size(ugibb,4) == 5
	[nx ny nz nu] = size(ugibb);
	if nu ~= 5, error size3d, end

	ugibb(:,:,:,[3 4]) = ugibb(:,:,:,[3 4]) * sqrt(2);

	total = squeeze(sum(ugibb,4));

	zx = zeros(nx,1);
	zy = zeros(1,ny);

	tmp = zeros(nx,ny,nz);
	tmp(1:nx-1,:,:) = ugibb(2:nx,:,:,1);
%	total = total + [ugibb(2:nx,:,1); zy];
	total = total + tmp;

	tmp = zeros(nx,ny,nz);
	tmp(:,1:ny-1,:) = ugibb(:,2:ny,:,2);
%	total = total + [ugibb(:,2:ny,2), zx];
	total = total + tmp;

	tmp = zeros(nx,ny,nz);
	tmp(2:nx,1:ny-1,:) = ugibb(1:nx-1,2:ny,:,3);
	total = total + tmp;

	tmp = zeros(nx,ny,nz);
	tmp(1:nx-1,1:ny-1,:) = ugibb(2:nx,2:ny,:,4);
	total = total + tmp;

	tmp = zeros(nx,ny,nz);
	tmp(:,:,1:nz-1) = ugibb(:,:,2:nz,5);
	total = total + tmp;
end

%
% 3D case - old way
%
if ndims(ugibb) == 4 && size(ugibb,3) == 5
	printm 'warn: old ugibb ordering'
	[nx ny nu nz] = size(ugibb);
	if nu ~= 5, error size3d, end

	ugibb(:,:,[3 4],:) = ugibb(:,:,[3 4],:) * sqrt(2);

	total = squeeze(sum(ugibb,3));

	zx = zeros(nx,1);
	zy = zeros(1,ny);

	tmp = zeros(nx,ny,nz);
	tmp(1:nx-1,:,:) = ugibb(2:nx,:,1,:);
%	total = total + [ugibb(2:nx,:,1); zy];
	total = total + tmp;

	tmp = zeros(nx,ny,nz);
	tmp(:,1:ny-1,:) = ugibb(:,2:ny,2,:);
%	total = total + [ugibb(:,2:ny,2), zx];
	total = total + tmp;

	tmp = zeros(nx,ny,nz);
	tmp(2:nx,1:ny-1,:) = ugibb(1:nx-1,2:ny,3,:);
	total = total + tmp;

	tmp = zeros(nx,ny,nz);
	tmp(1:nx-1,1:ny-1,:) = ugibb(2:nx,2:ny,4,:);
	total = total + tmp;

	tmp = zeros(nx,ny,nz);
	tmp(:,:,1:nz-1) = ugibb(:,:,5,2:nz);
	total = total + tmp;
end

%
% 2d case
%
if ndims(ugibb) == 3
	[nx,ny,nz] = size(ugibb);
	if nz ~= 4, error size2d, end

	ugibb(:,:,[3 4]) = ugibb(:,:,[3 4]) * sqrt(2);

	total = sum(ugibb,3);

	zx = zeros(nx,1);
	zy = zeros(1,ny);

	total = total + [ugibb(2:nx,:,1); zy];
	total = total + [ugibb(:,2:ny,2), zx];

	tmp4 = zeros(nx,ny);
	tmp4(1:nx-1,1:ny-1) = ugibb(2:nx,2:ny,4);
	total = total + tmp4;

	tmp3 = zeros(nx,ny);
	tmp3(2:nx,1:ny-1) = ugibb(1:nx-1,2:ny,3);
	total = total + tmp3;
end
