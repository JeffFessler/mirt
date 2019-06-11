 function output = ir_apply_tridiag_inv(sub, diags, sup, rhs)
%function output = ir_apply_tridiag_inv(sub, diags, sup, rhs)
%|
%| Solve system of equations with tridiagonal matrix (for many RHS),
%| i.e., apply the inverse of a tridiagonal matrix to rhs.
%|
%| in
%|	sub	[n-1]	subdiagonal
%|	diags	[n]	diagonal entries
%|	sup	[n-1]	superdiagonal
%|	rhs	[n m]	right hand side(s)
%|
%| out
%|	output	[n m]	T \ rhs
%|
%| 2013-09-15 Mai Le, University of Michigan
%| 2013-09-16 test, comments, white space added by JF
%| 2015-08-13 multiple RHS
%| 2016-01-05 timing test

if nargin == 1 && streq(sub, 'test'), ir_apply_tridiag_inv_test, return, end
if nargin == 1 && streq(sub, 'time'), ir_apply_tridiag_inv_time, return, end
if nargin < 4, ir_usage, end

n = length(diags);
m = size(rhs,2);
assert((length(sub) == n-1) & (length(sup) == n-1) & (size(rhs,1) == n), ...
	'vector lengths incompatible');

new_sup = zeros(n, 1, class(rhs));
new_arg = zeros(n, m, class(rhs));

new_sup(1) = sup(1) / diags(1);
new_arg(1,:) = rhs(1,:) / diags(1);
for ii = 2:n-1
	new_sup(ii) = sup(ii) / (diags(ii) - new_sup(ii-1) * sub(ii-1));
	new_arg(ii,:) = (rhs(ii,:) - new_arg(ii-1,:) * sub(ii-1)) ...
		/ (diags(ii) - new_sup(ii-1) * sub(ii-1));
end
new_arg(n,:) = (rhs(n,:) - new_arg(n-1,:) * sub(n-1)) ...
	/ (diags(n) - new_sup(n-1) * sub(n-1));

output = zeros(n, m, class(rhs));
output(n,:) = new_arg(n,:);
for ii = n-1:-1:1
	output(ii,:) = new_arg(ii,:) - new_sup(ii) * output(ii+1,:);
end


% ir_apply_tridiag_inv_test
% compare to \ as a simple test
function ir_apply_tridiag_inv_test

rng(7)
n = 18;
d = 3 + randn(n, 1);
sub = randn(n-1, 1);
%sup = randn(n-1, 1);
sup = sub; % make symmetric - usual case of interest for us
tri = diag(d) + diag(sub,-1) + diag(sup,1);
tri = single(tri);
%eig(tri)
if im
	pr cond(single(tri)) % Inf in octave @3.8.2_13+atlas+gcc5+x11
	pr cond(double(tri)) % 7.12
end

m = 16; % test with multiple RHS
rhs = randn(n,m) + 1i * rand(n,m); % stress with complex data
rhs = single(rhs);
xa = ir_apply_tridiag_inv(sub, d, sup, rhs);
xc = double(tri) \ double(rhs);
equivs(xa, xc)

if ~ir_is_octave % 3.8.2 \ has issues with single precision
	xb = tri \ rhs;
	try
	equivs(xa, xb)
	catch
	keyboard
	end
end

%jf_equal(xa, xb)


% ir_apply_tridiag_inv_time
% race vs \ and fft
function ir_apply_tridiag_inv_time

n = 2^9;
e2 = ones(n,1);
d2 = 3 * ones(n,1);
%Ts = spdiags([sub d sup], -1:1, n, n);
Ts2 = spdiags(-e2, -1, n, n) + spdiags(-e2, 1, n, n) + spdiags(d2, 0, n, n);
Tf2 = full(Ts2);
%Ts1 = single(Ts2); % not supported by matlab
Tf1 = single(Tf2);

sub = single(-e2(1:n-1));
sup = single(-e2(1:n-1));
d1 = single(d2);

om = [0:n-1]'/n*2*pi;
Hi = single(1 ./ (d2(1) - 2 * cos(om)));

m = n+1; % test with multiple RHS
rng(7);
rhs2 = randn(n,m) + 1i * rand(n,m); % stress with complex data
rhs1 = single(rhs2);
Hi_rep = repmat(Hi, [1 m]);

if 0 && im
	pr cond(single(Tf1))
	pr cond(double(Tf2))
	whos
end

nrep = 4; % repeat multiple times for "warm up"

for ii=1:nrep
	cpu etic
	xa = ir_apply_tridiag_inv(sub, d1, sup, rhs1);
	cpu etoc apply1
end

for ii=1:nrep
	cpu etic
	xf = ifft(Hi_rep .* fft(rhs1));
	cpu etoc fft1__
end

if 0 % much slower in double!
	for ii=1:nrep
		cpu etic
		xa = ir_apply_tridiag_inv(sub, d1, sup, rhs2);
		cpu etoc apply2
	end
end

for ii=1:nrep
	cpu etic
	x2s = Ts2 \ rhs2;
	cpu etoc back2s
end

for ii=1:nrep
	cpu etic
	x2f = Tf2 \ rhs2;
	cpu etoc back2f
end

equivs(xa, x2s)
equivs(xa, x2f)

%equivs(xa, xf) % not same due to periodic conditions
% clf, plot(real([xa(:,1)-xf(:,1)]), '.')
if 0
	im plc 1 2
	im(1, xa)
	im(2, xf)
end
