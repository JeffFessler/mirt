 function knot = ir_bspline_knot(ti, k)
%function knot = ir_bspline_knot(ti, k)
%| Given a set of sample point locations "ti",
%| determine reasonable knots for B-spline of order k.
%| Useful especially for zero boundary conditions?
%|
%| in
%|	ti [N]		sample locations
%|	k		order (nonnegative integer)
%|
%| out
%|	ti [N+k+1]	knot locations
%|
%| 2013-09-08, Jeff Fessler, University of Michigan

if nargin == 1 && streq(ti, 'test'), ir_bspline_knot_test, return, end
if nargin < 2, ir_usage, end

siz = size(ti);
ti = ti(:)'; % [1 N]

if rem(k,2) % odd
	d1 = ti(2) - ti(1);
	dN = ti(end) - ti(end-1);
	tmp = 1:((k+1)/2);
	knot = [ti(1) - d1*fliplr(tmp), ti, ti(end) + dN*tmp];
else % even
	mid = (ti(1:end-1) + ti(2:end)) / 2; % mid points
	d1 = mid(1) - ti(1);
	dN = ti(end) - mid(end);
	tmp = 1:(k/2+1);
	knot = [ti(1) - d1*fliplr(tmp), mid, ti(end) + dN*tmp];
end

if siz(1) ~= 1
	knot = knot';
end


function ir_bspline_knot_test
ti = [2 3 5 8 10]; % sample locations
t = linspace(-1, 14, 101);
t = [t, ti+1e-3, ti-1e-3];
t = sort(t);

im clf
for k=0:5
	knot = ir_bspline_knot(ti, k);
	out = ir_bspline_basis(t, k, knot);
	if im
		subplot(611 + k)
		plot(t, out, '-', t, sum(out,2), ':', ...
			ti, 0*ti, 'bo', knot, 0*knot, 'r.')
		axis([minmax(t)' 0 1.2])
		ytick([0 1])
		xtick(ti)
		ylabelf('k = %d', k)
	end
end
