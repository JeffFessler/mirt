 function tk = ir_chebyshev_poly(n)
%function tk = ir_chebyshev_poly(n)
%|
%| based on ChebyshevPoly.m by David Terr, Raytheon, 5-10-04
%|
%| Given nonnegative integer n, compute the
%| Chebyshev polynomial T_n. Return the result as a column vector whose mth
%| element is the coefficient of x^(n+1-m).
%| polyval(ir_chebyshev_poly(n), x) evaluates T_n(x).

if nargin < 1, ir_usage, end
if streq(n, 'test'), ir_chebyshev_poly_test, return, end

if n==0
	tk = 1;
elseif n==1
	tk = [1 0]';
else
	tkm2 = zeros(n+1,1);
	tkm2(n+1) = 1;
	tkm1 = zeros(n+1,1);
	tkm1(n) = 1;

	for k=2:n
		tk = zeros(n+1,1);

		for e=n-k+1:2:n
			tk(e) = 2*tkm1(e+1) - tkm2(e);
		end

		if mod(k,2)==0
			tk(n+1) = (-1)^(k/2);
		end

		if k<n
			tkm2 = tkm1;
			tkm1 = tk;
		end
	end
end


% ir_chebyshev_poly_test()
function ir_chebyshev_poly_test
nlist = [0:4];
x = 1 * linspace(-1,1,101)';
y = zeros(numel(x), numel(nlist));
for in=1:numel(nlist);
	n = nlist(in);
	p = ir_chebyshev_poly(n);
	y(:,in) = polyval(p,x);
end
if im
	clf
	plot(x, y)
	axis([-1 1 -1.1 1.1])
	xtick([-1 0 1]), ytick([-1 0 1])
	titlef('Chebychev Polynomials')
end
