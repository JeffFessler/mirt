 function [basis, Hw] = deconv_z(psf, fun, t)
% psf: 3-tap symmetric psf: [a _b_ a]
% H(z) = 1 / (az + b + a/z)

if numel(psf) == 1
	basis = fun(t) / psf;
	Hw = @(t) psf * ones(size(t));
return
end

if length(psf) ~= 3, error 'not done', end
a = psf(1);
b = psf(2);
c = b / a / 2;
p = -c + sign(c) * sqrt(c^2 - 1); % pole

scale = 1/a * 1/(p - 1/p);
basis = fun(t);
for n=1:9
	basis = basis + p^n * (fun(t-n) + fun(t+n));
end
basis = scale * basis;

Hw = @(om) 1 ./ (b + 2 * a * cos(om));
