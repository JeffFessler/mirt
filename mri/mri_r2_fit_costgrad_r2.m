  function [cost grad] = mri_r2_fit_costgrad_r2(r2, data)
%|function [cost grad] = mri_r2_fit_costgrad_r2(r2, data)
%| for R2=1/T2 estimation from images

if nargin == 1 && streq(r2, 'test'), mri_r2_fit_costgrad_r2_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

yi = data.yi; % [(nd) nt] images
pd = data.pd; % [(nd)] pd maps
te = data.te; % echo times

[yb yi] = mri_r2_fit_mean(pd, r2, yi, te);
cost = sum(abs(yi(:) - yb(:)).^2) / 2;

grad = mri_r2_fit_grad_r2(pd, r2, yi - yb, te);

if isfield(data, 'R')
	R = data.R;
	cost = cost + R.penal(R, r2);
	grad = grad + R.cgrad(R, r2);
end

grad = reshape(grad, size(r2));


% mri_r2_fit_mean()
function [yb yi] = mri_r2_fit_mean(pd, r2, yi, te)

nt = numel(te);
yi = reshapee(yi, [], nt); % [*nd nt]
yb = zeros(size(yi)); % [*nd nt]
for it=1:nt
	yb(:,it) = pd(:) .* exp(-te(it) * r2(:));
end


% mri_r2_fit_grad_r2()
function grad = mri_r2_fit_grad_r2(pd, r2, resid, te)

nt = numel(te);

grad = 0;
for it=1:nt
	tmp = pd(:) .* exp(-te(it) * r2(:));
	grad = grad + te(it) * real(conj(tmp) .* resid(:,it));
end


function mri_r2_fit_costgrad_r2_test

r2 = 0.01;
te = [0:3] * 20;
pd = 2;
yb = pd .* exp(-te * r2);
rng(0)
yi = yb + 0.1 * randn(size(yb));

%R = Reg1(true(1), 'beta', 0, 'order', 0, 'type_penal', 'mat');

data.yi = yi;
data.pd = pd;
data.te = te;
%data.R = R;

r2list = linspace(0, 2*r2, 101);
cost = zeros(size(r2list));
for ir=1:numel(r2list)
	cost(ir) = mri_r2_fit_costgrad_r2(r2list(ir), data);
end

r0 = 0.004;
[cost0 grad0] = mri_r2_fit_costgrad_r2(r0, data);

if im
	clf
	plot(r2list, cost, '-', r2list, cost0 + grad0*(r2list-r0), '--')
end
