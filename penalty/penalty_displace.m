 function displace = penalty_displace(offsets, sizes)
%function displace = penalty_displace(offsets, sizes)
%|
%| Convert scalar offsets to vector displacements, i.e.,
%| find d1,d2,... such that d1*1 + d2*n1 + d3*n1*n2 + ... = offset.
%| Example: if offset = n1, then [d1 d2 d3] = [0 1 0]
%| Assumes that d <= floor(n/2) for n > 2
%|
%| in
%|	offsets	[LL]		see penalty_offsets.m
%|	sizes	[1 ndim]
%|
%| out
%|	displace [LL ndim]	[d1 d2 ...] for each offset
%|
%| Copyright 2006-12-6, Jeff Fessler, University of Michigan

if nargin == 1 && streq(offsets, 'test'), penalty_displace_test, return, end
if nargin < 2, ir_usage, end

ndim = numel(sizes);
offsets = offsets(:); % [LL 1]
displace = nan(length(offsets), ndim);


if 1 % simple way that assumes d / n < 1/2, which is fine for n > 2
	if length(sizes) == 2 && sizes(1) == 2 % [2 ny] special case!
		if isequal(offsets, [1 2 3 1]') % [1 nx nx+1 nx-1]
			displace = [1 0; 0 1; 1 1; -1 1];
		elseif isequal(offsets, [1]) % [1]
			displace = [1 0];
			warn 'caution: offset=1 is ambiguous when n1=2'
		elseif isequal(offsets, [2]) % [n1]
			displace = [0 1];
		elseif isequal(offsets, [1 2]') % [2 n2]
			displace = [1 0; 0 1];
		else
			pr offsets
			fail 'not done'
		end

	elseif any(sizes(1:end-1) == 2) 
		fail 'not done: ambiguous'

	else
		residual = offsets;
		for id=ndim:-1:1
			nd = prod(sizes(1:id-1));
			displace(:,id) = round(residual / nd);
			residual = residual - displace(:,id) * nd;
		end
		if any(residual ~= 0), fail 'bug', end
	end

else

	half = max(floor(sizes/2), 1); % upper bound on dx is almost half the size
	%half = sizes/2;

	for ll=1:length(offsets)
		offset = double(offsets(ll)); % trick: necessary!
		subval = 0;
		dis = zeros(1,ndim);
		for id=ndim:-1:2
			tmp = offset + sum((half(1:id-1)-1) .* [1 sizes(1:id-2)]);
			tmp = tmp + prod(sizes(1:id)) - subval;
			dis(id) = floor(tmp / prod(sizes(1:id-1))) - sizes(id);
			subval = subval + dis(id) * prod(sizes(1:id-1));
		end
		dis(1) = offset - subval;
	
%		if ~(all(abs(dis) < half || (dis == 1 && half == 1)))
		if ~all(abs(dis) <= half)
			pr dis
			pr half
			warn 'bug?'
		end
		displace(ll,:) = dis;
	end
end

offset_check = displace * [1 cumprod(sizes(1:end-1))]';
jf_equal(offsets, offset_check)

end % penalty_displace


% penalty_displace_test()
function penalty_displace_test

jf_equal(penalty_displace([], [1]), nan(0,1))
jf_equal(penalty_displace([], [1 1]), nan(0,2))

if 1 % check tiny
	% when nx=ny=2 then [1 nx nx+1 nx-1] = [1 2 3 1] so it is ambiguous!
	nx = 2; ny = 2;
	offsets = [1 nx nx+1 nx-1];
	dd = penalty_displace(offsets, [nx ny]);
	jf_equal(dd, [1 0; 0 1; 1 1; -1 1])

	nx = 2; ny = 3;
	dd = penalty_displace(offsets, [nx ny]);
	jf_equal(dd, [1 0; 0 1; 1 1; -1 1])
end

nx = 3; ny = 3;
offsets = [1 nx nx+1 nx-1];
dd = penalty_displace(offsets, [nx ny]);
jf_equal(dd, [1 0; 0 1; 1 1; -1 1])

nx = 100; ny = 80; % 2d
[ix iy] = ndgrid(-2:2, -2:2);
offsets = col(ix + iy * nx);
dd = penalty_displace(offsets, [nx ny]);
jf_equal(dd, [ix(:) iy(:)])

nx = 10; ny = 8; nz = 7; % 3d
[ix iy iz] = ndgrid(-2:2, -2:2, -2:2);
offsets = col(ix + iy * nx + iz * nx * ny);
dd = penalty_displace(offsets, [nx ny nz]);
jf_equal(dd, [ix(:) iy(:) iz(:)])

dim = [8 9]; % for won
offsets2 = penalty_offsets('2d:hvd', dim);
dd2 = penalty_displace(offsets2, dim);

for nz=1:3
	dim3 = [dim nz];
	dd3 = penalty_displace(offsets2, dim3);
	jf_equal(dd3, [dd2 zeros(4,1)])
end

if 1
	dim3 = [10 20 30];
	offsets = penalty_offsets('3d:26', dim3);
	[[1:13]', penalty_displace(offsets, dim3)]
end

end % penalty_displace_test()
