  function [grad, curv, fim, diml] = de_wls_dercurv(data, sil, curvtype, iblock, nblock)
%|function [grad, curv, fim, diml] = de_wls_dercurv(data, sil, curvtype, iblock, nblock)
%|
%| evaluate gradients and curvatures for DE WLS data-fit function:
%| U(x) = \sum_i U_i(s)
%| U_i(s) = \sum_m=1^M U_im(s)
%| U_im(s_i) = 1/2 w_im (h_im - f_m(s_i))^2
%|
%| where f : R^L -> R^M is associated with polyenergetic beam-hardening
%| for a single material type (usually water).
%|
%| in
%|	data	{him, wim, ftab}
%|				him is [nb na M]
%|				wim is same, or just "1"
%|	sil	[nb na L]	material density line integrals
%|	curvtype		todo
%|	iblock, nblock		for OS type methods
%|
%| out
%|	grad	[nb na L]	\dot hi(l)
%|	curv	[nb na L]	L-separable surrogate curvature for each U_i
%|	fim	[nb na M]	f_m(s_i), should approximate yi
%| ?	dil	[nb na ?]	\dot f_m(s_i)
%|
%| todo: this could generalize to any nonlinear functions g!
%|
%| Copyright 2008-9-27, Jeff Fessler

if nargin < 2, help(mfilename), error args, end

[him wim ftab] = deal(data{:});

fit = ftab.fit;

size_sil = size(sil);
%if fit.LL == 1
%	if ns(end) == 1, ns = ns(1:end-1); end
%else
%	ns = ns(1:end-1);
%end

if 1 == size(sil,2)
	sil = reshapee(sil, [], fit.LL);
end

if 1 == size(him,2)
	him = reshapee(him, [], fit.MM);
end

if nargin == 5 % handle subset case
	ia = iblock:nblock:size(him,2);
	him = him(:,ia,:);
	if length(wim) > 1
		wim = wim(:,ia,:);
	end
	warn 'not tested'
end

%printm fim
fim = fit.fmfun(sil);
%fim = reshape(fim, size(him));

%printm ci
%ci = fit.ls_curv(sil, fim, him); % [() 1]
%ci = reshape(ci, size(him)); % todo?

%printm di
%dilm = fit.fgrad(sil); % [() L M]
%di = reshape(di, size(him));

[grad curv] = de_ftab_wls_grad_curv(fit, him, fim, sil, wim); % [*ns L]
grad = reshape(grad, size_sil);
curv = reshape(curv, size_sil);

if any(isnan(grad(:))) || any(isnan(curv(:)))
	warn 'nan'
	keyboard
end

%deriv = fit.ls_grad(yim, sil); % [() L M]
%deriv = wi .* di .* (fi - yi);
%curv = wi .* ci;

%
% de_ftab_wls_grad_curv()
% gradient of WLS cost function
% in
%	fit	strum
%	him	[(ns) M]	\hat{f}_im (log estimates)
%	fim	[(ns) M]	f_im (log predictions from sil), possibly empty
%					fit.fmfun(sil)
%	sil	[(ns) L]	component material density integrals
%	wim	[(ns) M]	weights
% out
%	grad	[(ns) L]	gradient of WLS cost function
%	curv	[(ns) L]	curvatures of WLS cost function
%
function [grad curv] = de_ftab_wls_grad_curv(fit, him, fim, sil, wim)

M = fit.MM;
L = fit.LL;
ns = size(sil);
if L == 1
	if ns(end) == 1, ns = ns(1:end-1); end
else
	if fit.LL ~= ns(end), error 'sil size', end; ns = ns(1:end-1);
end

sil = reshapee(sil, [], L); % [*ns L]
him = reshapee(him, [], M); % [*ns M]

if isempty(fim)
	fim = fit.fmfun(sil); % [*ns M]
else
	fim = reshapee(fim, [], M); % [*ns M]
end

fgrad = fit.fgrad(sil); % [*ns L M]

if ~isvar('wim') || isempty(wim)
	wim = ones(size(fim));
else
	wim = reshapee(wim, [], M); % [*ns M]
end

if ~streq(fit.ctype, 'pre10')
	fail 'not done'
end

tmp = fit.mac_eff; % [M L]
curv = abs(tmp)' * (abs(tmp) * ones(L,1)); % [L 1]
curv = max(wim,[],2) * curv'; % [*ns L]

grad = 0;

for mm=1:M
	err = fim(:,mm) - him(:,mm); % [*ns 1]
	tmp = repmat(wim(:,mm) .* err, [1 L]); % [*ns L]
	grad = grad + tmp .* fgrad(:,:,mm);
end

grad = reshape(grad, [ns L]); % [(ns) L]
curv = reshape(curv, [ns L]); % [(ns) L]
