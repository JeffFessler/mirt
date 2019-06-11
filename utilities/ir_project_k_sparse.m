 function y = ir_project_k_sparse(x, k)
%function y = ir_project_k_sparse(x, k)
%|
%| zero out all but k largest |entries| in each column of x
%|
%| in
%|	x	[n m]
%|	k	[1]		natural number
%|
%| out
%|	x	[n m]
%|
%| 2013-12-10 Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), ir_project_k_sparse_test, return, end
if nargin < 2, ir_usage, end

[s i1] = sort(abs(x), 1, 'descend');
i1 = i1(1:k,:);
i2 = 1:size(x,2);
i2 = repmat(i2, [k 1]);
sub = sub2ind(size(x), i1, i2);

y = zeros(size(x), class(x));
y(sub) = x(sub);

function ir_project_k_sparse_test
x = reshape(1:12, 3, 4);
x = [magic(4) magic(4)];
k = 3;
y = ir_project_k_sparse(x, k);
if im
	pr x
	pr y
end
