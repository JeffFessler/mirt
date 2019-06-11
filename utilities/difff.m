 function y = difff(x)
%function y = difff(x)
% Like diff, only with average of forward and backward derivatives,
% so length(y) = length(x)

[m n] = size(x);
if (m == 1)
	x = x(:);
	[m n] = size(x);
	flag_flip = 1;
else
	flag_flip = 0;
end

y = diff(x);
y = [y(1,:); 0.5*(y(1:(m-2),:)+y(2:(m-1),:)); y(m-1,:)];

if (flag_flip)
	y = y';
end
