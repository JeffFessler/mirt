 function out = jf_mid3(pn, zz, dim)
%function out = jf_mid3(pn, zz, dim)
%|
%| Extract "middle" slices in all 3 dimensions of 3d array 'zz'
%| or the middle slice from dimension 'dim'.
%| For an even-sized dimension, it uses slice n/2+1 (per FFT), as of 2022-08-31.

if nargin < 2, ir_usage, end

if streq(zz, 'test'), jf_mid3_test(pn); out = []; return, end

assert(ndims(zz) == 3)
[nx, ny, nz] = size(zz);

%jf_mid3_fun = @(n) ceil(n/2); % before 2022-08-31
jf_mid3_fun = @(n) ceil((n+1)/2);
i1 = jf_mid3_fun(nx);
i2 = jf_mid3_fun(ny);
i3 = jf_mid3_fun(nz);

if nargin > 2 && ~isempty(dim)
	nn = size(zz, dim);
	switch dim
	case 1
		out = zz(i1,:,:);
	case 2
		out = zz(:,i2,:);
	case 3
		out = zz(:,:,i3);
	otherwise
		fail('bad dim %d', dim)
	end
return
end

xy = zz(:,:,i3);
xz = permute(zz(:,i2,:), [1 3 2]);
zy = permute(zz(i1,:,:), [2 3 1])';
out = [xy, xz; zy, zeros(nz,nz)];
end


function jf_mid3_test(pn)
	x = rand(3,4,5);
	y = jf_mid3(pn, x);
	assert(all(size(y) == [8 9]))
	y = jf_mid3(pn, x, 1);
	assert(all(y == x(2,:,:), 'all'))
	y = jf_mid3(pn, x, 2);
	assert(all(y == x(:,3,:), 'all'))
	y = jf_mid3(pn, x, 3);
	assert(all(y == x(:,:,3), 'all'))
end
