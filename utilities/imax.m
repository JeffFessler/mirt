 function [ii, i2] = imax(a, flag)
%function [ii, i2] = imax(a, flag)
% Return index of maximum of each column of a
% flag = 2 to treat as 2d
% flag = 3 to treat as 3d

if nargin < 1, ir_usage, end

if nargin == 1
	[dum, ii] = max(a);
	return
end

%
% multidimensional cases
%

if flag == 2
	if ndims(a) ~= 2
		error 'not done'
	else
		[dum, ii] = max(a(:));
		[i1, i2] = ind2sub(size(a), ii);
		if nargout <= 1
			ii = [i1 i2];
		elseif nargout == 2
			ii = i1;
		else
			error '1 or 2 outputs only'
		end
	end

elseif flag == 3
	if ndims(a) ~= 3
		error 'not done'
	else
		[dum, ii] = max(a(:));
		[i1, i2, i3] = ind2sub(size(a), ii);
		ii = [i1 i2 i3];
	end
end
