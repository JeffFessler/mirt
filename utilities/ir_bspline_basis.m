 function out = ir_bspline_basis(t, k, ti)
%function out = ir_bspline_basis(t, k, ti)
%|
%| Return B-spline basis of order k >= 0, for knots ti, at locations t.
%| This is useful for didactic plotting and small problems
%| but is not efficient (compute or memory) for large problems.
%|
%| in
%|	t [M 1]		values at which to evaluate B-spline
%|	k		order (nonnegative integer)
%|	ti [N+k+1]	knot locations
%|
%| out
%|	out [M N]	B_{i,k}(t) for each t and for i=1,...,N
%|
%| 2013-09-08, Jeff Fessler, University of Michigan

if nargin == 1 && streq(t, 'test'), ir_bspline_basis_test, return, end
if nargin < 3, ir_usage, end

t = t(:); % [M 1]
ti = ti(:)'; % [1 N]
siz = size(t);

N = numel(ti) - k - 1;
M = numel(t);

if any(diff(ti) <= 0)
	fail 'knot locations must be strictly increasing'
end

if k == 0
	tmp = t(:); % [*M 1]
	left = ti(1:N);
	right = ti(2:N+1);
	out = (outer_sum(tmp, -left) >= 0) & (outer_sum(-tmp, right) > 0);
	out = single(out);
else
	b1 = ir_bspline_basis(t, k-1, ti(1:end-1));
	b2 = ir_bspline_basis(t, k-1, ti(2:end));
	out = zeros(M, N, 'single');
	for ii=1:N
		s1 = (t - ti(ii)) / (ti(ii+k) - ti(ii));
		s2 = (ti(ii+k+1) - t) / (ti(ii+k+1) - ti(ii+1));
		f1 = s1 .* b1(:,ii);
		f2 = s2 .* b2(:,ii);
		out(:,ii) = 1 * f1 + 1 * f2;
	end
end


function ir_bspline_basis_test
k = 3;
ti = [2 3 5 8 9 11 12 14 15];
%ti = [2:12];
t = linspace(0, 17, 101);
t = [t, ti+1e-3, ti-1e-3];
t = sort(t);

im clf
for k=0:3
	out = ir_bspline_basis(t, k, ti);
	subplot(411 + k)
	plot(t, out, '-', t, sum(out,2), ':', ti, 0*ti, 'bo')
	axis([minmax(t)' 0 1.2])
	ytick([0 1])
end
