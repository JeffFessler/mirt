function out = jf_mid3(pn, zz, dim)
% extract middle slices in all 3 dimensions of 3d array 'zz'
% or the middle slice from dimension 'dim'
if nargin < 2, ir_usage, end

if nargin > 2 && ~isempty(dim)
	nn = size(zz, dim);
	switch dim
	case 1
		out = zz(ceil(end/2),:,:);
	case 2
		out = zz(:,ceil(end/2),:);
	case 3
		out = zz(:,:,ceil(end/2));
	otherwise
		fail('bad dim %d', dim)
	end
return
end

nz = size(zz,3);
xy = zz(:,:,ceil(end/2));
xz = permute(zz(:,ceil(end/2),:), [1 3 2]);
zy = permute(zz(ceil(end/2),:,:), [2 3 1])';
out = [ xy, xz; zy, zeros(nz,nz)];
