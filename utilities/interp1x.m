 function yi = interp1x(x, y, xi, flag1)
%function yi = interp1x(x, y, xi, flag1)
% Given ordered x,y pairs and values xi, interpolate.
% if flag1 is used, then use linear interpolation outside range.
% by default, only 0 order interpolation is used.

if 0
	[x,j] = sort(x(:));
	y = y(j);
	t = diff(x) == 0;
	if any(t)
		t = find(t)
		x(t) = [];
		y(t) = [];
	end
	yi = 0*xi;
	t = xi < x(1);
	if sum(t),	yi(t) = x(1) * ones(sum(t),1);	end
	t = xi > x(length(x));
	if sum(t),	yi(t) = x(length(x)) * ones(sum(t),1);	end
	t = xi >= x(1) & xi <= x(length(x));
	yi(t) = interp1(x, y, xi(t));
%	yi(t) = table1([x y], xi(t));
return
end

if (nargin < 3), help interp1x, return, end

N = length(x);
M = numel(xi);
yi = xi;	x=x(:); y=y(:); xi=xi(:);

if (N ~= length(y)), error(['ERROR: x,y must be same length']), end

[x,i] = sort(x);
y = y(i);

ii = sum(x*ones(1,M) < ones(N,1)*xi', 1);
ig = ii>=1 & ii < N;
if nargin == 4
	yi(ii==0) = y(1) + (xi(ii==0) - x(1))/(x(2)-x(1))*(y(2)-y(1));
	yi(ii>=N) = y(N-1) + (xi(ii>=N) - x(N-1))/(x(N)-x(N-1))*(y(N)-y(N-1));
else
	yi(ii==0) = y(1)*ones(sum(ii==0),1);
	yi(ii>=N) = y(N)*ones(sum(ii>=N),1);
end

ii = ii(ig);
yi(ig) = y(ii) + (xi(ig)-x(ii))./(x(ii+1)-x(ii)) .* (y(ii+1) - y(ii));
