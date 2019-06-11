  function [ii i2] = imin(a, flag2d)
%|function [ii i2] = imin(a, flag2d)
%| Return index of minimum of each column of a
%| see imax()

if nargin < 1, ir_usage, end

if nargin == 1
	[dum ii] = min(a);
else
	if nargout <= 1
		ii = imax(-a, flag2d);
	elseif nargout == 2
		[ii i2] = imax(-a, flag2d);
	else
		fail 'too many outputs'
	end
end
		
