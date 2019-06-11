 function data = poisson1(xmean, seed)
%function data = poisson1(xmean, seed)
%| Generate poisson random column vector with mean xmean
%| by summing exponentials.
%| Optional seed for rng().
%| This is efficient only for small mean values, eg < 20.

if nargin < 1, ir_usage, end
if streq(xmean, 'test'), poisson1_test, return, end

if isvar('seed') && ~isempty(seed)
	rng(seed)
end

dim = size(xmean);

data = zeros(dim);
i_do = ones(dim);
ee = exp(xmean);

while any(i_do(:))
	i_do = ee >= 1;
	data(i_do) = 1 + data(i_do);
	ee = ee .* rand(dim) .* i_do;
end

data = data - 1;


% poisson1_test()
function poisson1_test
n = 2^8;
t = reshape(linspace(1, 30, n^2), n, n);
cpu etic
poisson1(t);
cpu etoc 'fessler poisson time'
if exist('poissrnd') == 2
	cpu etic
	poissrnd(t);
	cpu etoc 'matlab poissrnd time'
end

tmp = 9 * ones(n^2,1);
tmp = poisson1(tmp, 7);

if im
	pr '[mean(tmp) var(tmp)]'
end
