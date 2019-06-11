  function offset_zxy = reg_offset_xyz_to_zxy(offset_xyz, dim_xyz)
%|function offset_zxy = reg_offset_xyz_to_zxy(offset_xyz, dim_xyz)
%|
%| convert usual offsets for xyz image ordering
%| to offsets appropriate for zxy image ordering
%|
%| in
%|	offset_xyz	[1 noffset]
%|	dim_xyz		[3]
%|
%| out
%|	offset_zxy	[1 noffset]
%|
%| Copyright 2009-5-3, Jeff Fessler, University of Michigan

if nargin == 1 && streq(offset_xyz, 'test')
	reg_offset_xyz_to_zxy_test
return
end
if nargin < 2, help(mfilename), error(mfilename), end

if numel(dim_xyz) ~= 3
	fail('expected dim_zxy to be [3]')
end

if ischar(offset_xyz)
	offset_xyz = penalty_offsets(offset_xyz, dim_xyz);
end

dd_xyz = penalty_displace(offset_xyz, dim_xyz); % [n 3]
dd_zxy = dd_xyz(:, [3 1 2]); % [n 3]

% ensure last nonzero displacement is positive, so offset is positive.
% this should not affect directional beta_l values.
if 1
	for io=1:nrow(dd_zxy)
		if dd_zxy(io,3) < 0
			dd_zxy(io,:) = -dd_zxy(io,:);
		elseif dd_zxy(io,3) == 0 && dd_zxy(io,2) < 0
			dd_zxy(io,:) = -dd_zxy(io,:);
		end
	end
end

dim_xyz = dim_xyz(:)'; % [1 3]
dim_zxy = dim_xyz([3 1 2]); % [1 3]
offset_zxy = [1, cumprod(dim_zxy(1:end-1))] * dd_zxy'; % [1 n] = [1 3] * [3 n]

if any(offset_zxy <= 0)
	warn 'negative offsets?'
	pr offset_zxy
	keyboard
end


% reg_offset_xyz_to_zxy_test()
function reg_offset_xyz_to_zxy_test
dim_xyz = [10 20 30];
dim_zxy = dim_xyz([3 1 2]);
o_xyz = penalty_offsets('3d:26', dim_xyz);
o_zxy = penalty_offsets('3d:26', dim_zxy);
d_xyz = penalty_displace(o_xyz, dim_xyz)';
d_zxy = penalty_displace(o_zxy, dim_zxy)';
tmp = reg_offset_xyz_to_zxy(o_xyz, dim_xyz);
jf_equal(sort(tmp), sort(o_zxy))
