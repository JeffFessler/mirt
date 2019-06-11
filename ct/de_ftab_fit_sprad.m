 function sr = de_ftab_fit_sprad(fit, sl)
%function sr = de_ftab_fit_sprad(fit, sl)
% for each s vector compute spectra radius of hessian of LS fit
% in
%	sl	cell{L}	[(Nd),L]
% out
%	sr	[(Nd),M]
% applies only to 'exp' fit

sll = ndgrid_jf('mat', sl{:}); % [(Nd),L]
Nd = size(sll); LL = Nd(end); Nd = Nd(1:end-1);
sll = reshape(sll, [], LL); % [*Nd,L]

if LL ~= 2, error 'only L=2 done', end
MM = fit.MM;

%sr = cell(fit.MM,1);
sr = zeros(prod(Nd), MM);
for mm=1:MM
	alf = fit.coef{mm}; % [ne,1]
	mac = fit.mac{mm}; % [ne,L]
	g0 = mac' * alf; % gradient of fit at s=0 (largest point) 
	h0 = mac' * diag(alf) * mac - g0 * g0'
	norm(g0)
	norm(g0 * g0')
	norm(h0)
	for is=1:prod(Nd)
		ss = sll(is,:)'; % [L,1]
		q = alf .* exp(-mac * ss); % [ne,1]
		q = q / sum(q);
%		plot(q), drawnow
		g = mac' * q; % gradient of f_m(s)
		h = mac' * diag(q) * mac - g * g';
		sr(is,mm) = norm(g * g');
%		sr(is,mm) = norm(h);
	end
end
sr = reshape(sr, [Nd MM]);
