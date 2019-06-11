  function y = interp1_lagrange(xi, yi, x)
%|function y = interp1_lagrange(xi, yi, x)
%|
%| Lagrange (1D) interpolation: polynomial of order length(xi)
%|
%| in
%|	xi [n 1], yi [n 1], x [m 1]
%| out
%|	y [m 1]		if yi is empty, then "y" is the m,n interpolation matrix
%|
%| Copyright 2003-5-13, Jeff Fessler, University of Michigan

if nargin == 1 && streq(xi, 'test'), interp1_lagrange_test, return, end
if nargin < 3, ir_usage, end

n = length(xi);
%num = outer_sum(x, -xi);		% [m n] 
den = outer_sum(xi, -xi); 
den = den + (1 - den) .* eye(n);	% i \neq j
den = prod(den,2);

out = outer_sum(x, -xi);		% [m n] 

if isempty(yi)
	y = zeros(length(x), n);
	for ii=1:n
		t = out;
		t(:,ii) = [];			% [m n-1]
		y(:,ii) = prod(t,2) / den(ii);
	end

else
	yi = yi(:) ./ den;	% [n 1]

	y = zeros(size(x));
	for ii=1:n
		t = out;
		t(:,ii) = [];			% [m n-1]
		y = y + yi(ii) * prod(t, 2);
	end
end


function interp1_lagrange_test
xi = [1 3 4 7]';
%xi = [-2:2]';
yi = [6 2 5 1]';
%yi = xi == 0;
x = linspace(0,8,101)';
%x = linspace(-2,2,101)';
y1 = interp1_lagrange(xi, yi, x);
y2 = interp1_lagrange(xi, [], x);
if im
	plot(xi, yi, 'o', x, y1, 'y-', x, y2, ':', x, y2*yi, 'c.')
end
