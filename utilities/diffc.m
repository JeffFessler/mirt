  function y = diffc(x)
%|function y = diffc(x)
%|
%| circulant version of diff()

if nargin < 1, ir_usage, end

[m n] = size(x);
y = diff(x);
if (m == 1)
 	y(n) = x(1) - x(n);
else
	y(m,:) = x(1,:) - x(m,:);
end
