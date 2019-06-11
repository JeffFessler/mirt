  function [cost grad] = mri_r2_fit_costgrad_pd(pd, data)
%|function [cost grad] = mri_r2_fit_costgrad_pd(pd, data)
%| gradient wrt "proton density" (pd) image
%| for R2=1/T2 estimation from images

if nargin == 1 && streq(pd, 'test'), mri_r2_fit_costgrad_pd_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

yi = data.yi; % [(nd) nt] images
r2 = data.r2; % [(nd)] R2=1/T2 maps
te = data.te; % echo times

[yb yi] = mri_r2_fit_mean(pd, r2, yi, te);
cost = sum(abs(yi(:) - yb(:)).^2) / 2;

grad = mri_r2_fit_grad_pd(pd, r2, yi - yb, te);

if isfield(data, 'R')
	R = data.R;
	cost = cost + R.penal(R, pd);
	grad = grad + R.cgrad(R, pd);
end
grad = reshape(grad, size(pd));


% mri_r2_fit_mean()
function [yb yi] = mri_r2_fit_mean(pd, r2, yi, te)

nt = numel(te);
yi = reshapee(yi, [], nt); % [*nd nt]
yb = zeros(size(yi)); % [*nd nt]
for it=1:nt
	yb(:,it) = pd(:) .* exp(-te(it) * r2(:));
end


% mri_r2_fit_grad_pd()
function grad = mri_r2_fit_grad_pd(pd, r2, resid, te)

nt = numel(te);

grad = 0;
for it=1:nt
	tmp = exp(-te(it) * r2(:));
	grad = grad - tmp .* resid(:,it);
end


function mri_r2_fit_costgrad_pd_test

r2 = 0.01;
te = [0:3] * 20;
pd = 2;
yb = pd .* exp(-te * r2);
rng(0)
yi = yb + 0.1 * randn(size(yb));

%R = Reg1(true(1), 'beta', 3, 'order', 0, 'type_penal', 'mat');

data.yi = yi;
data.r2 = r2;
data.te = te;
%data.R = R;

pdlist = linspace(0, 4, 101);
cost = zeros(size(pdlist));
for ii=1:numel(pdlist)
	cost(ii) = mri_r2_fit_costgrad_pd(pdlist(ii), data);
end

p0 = 1.5;
[cost0 grad0] = mri_r2_fit_costgrad_pd(p0, data);

if im
	clf
	plot(pdlist, cost, '-', pdlist, cost0 + grad0*(pdlist-p0), '--')
end
