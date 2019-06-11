  function fit = de_ftab_fit_exp3(xrs, mac, varargin)
%|function fit = de_ftab_fit_exp3(xrs, mac, [options])
%|
%| Fit to DE table fm(), suitable for subsequent interpolation / extrapolation.
%| Uses a special experimental exponential model:
% todo -log(sum_k p_k exp(-m_k . s))
% fit a 3-term exponential model that has the proper derivatives
% at 0, and at +/- infinity
%|
%| in
%|	xrs	strum		X-ray spectra; see xray_read_spectra.m
%|	mac	[Ne,L]		mass attenuation coefficients
%|
%| option
%|	'show'	1|0		plot?
%|
%| out
%|	fit	strum		strum object for fitted f_m(s_1,...,s_L)
%|
%|	methods:
%|	.fmfun(sll)		fm function evaluation, for stacked array sll
%|	.fgrad(sll)		fm gradient evaluation, for stacked array sll
%|	.show_sp(en, sp)	plot true spectrum vs fitted spectrum
%|	.show_fm(sl, fm)	mesh plot of fm and its fit
%|	.show_err(sl, fm)	mesh plot of fit error
%|	.mac_eff		effective mass atten coef based on fit
%|				(valid for 'exp' only)
%|
%| Copyright 2008-8-10, Jeff Fessler, University of Michigan

%if nargin == 1 && streq(xrs, 'test'), de_ftab_fit_exp3_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.show = false;
arg = vararg_pair(arg, varargin);

LL = ncol(mac);
if LL ~= 1, fail 'only 1 component done', end

MM = ncol(xrs.Ide);

fit.MM = MM;
fit.type = 'exp3';

massbar = (xrs.Ide' * mac) ./ sum(xrs.Ide)'; % [M,1]

%if ~isempty(wt), warn 'weighting ignored', end
% if isempty(wt), wt = num2cell(ones(MM,1)); end

fit.kev = cell(1,MM);
fit.mac = cell(1,MM);
fit.coef = cell(1,MM);
for mm=1:MM
	good = xrs.sp(:,mm) > 0;
	sp_m = xrs.sp(good,mm); % nonzero spectrum
	sp_m = sp_m / sum(sp_m);
	[mac_m ii] = sort(mac(good));
	if any(mac_m(2:end)) == mac_m(1), error 'non unique min', end
	if any(mac_m(1:end-1)) == mac_m(end), error 'non unique max', end
	kev_m = xrs.en(good);
	kev_m = kev_m(ii);
	q1 = sp_m(1);
	qN = sp_m(end);
	b1 = mac_m(1);
	bN = mac_m(end);
	mac_mid = (massbar(mm) - q1 * b1 - qN * bN) / (1 - q1 - qN);
	fit.mac{mm} = [b1 mac_mid bN]';
	fit.coef{mm} = [q1 1-q1-qN qN]';

	kev_mid = interp1(mac_m, kev_m, mac_mid);
	fit.kev{mm} = [kev_m(1) kev_mid kev_m(end)]';
end

meth = {'fmfun', @de_ftab_fit_exp_eval, '(sll)'; ...
	'fgrad', @de_ftab_fit_exp_grad, '(sll)'; ...
	'show_err', @de_ftab_fit_show_err, '(sl, fm)'; ...
	'show_fm', @de_ftab_fit_show_fm, '(sl, fm)'; ...
	'show_sp', @de_ftab_fit_show_sp, '(en, sp)'; ...
	};
fit = strum(fit, meth);

if arg.show
	fit.show_fm(sl, fm);
end

end % de_ftab_fit_exp3()


%
% de_ftab_fit_exp()
% todo: more documentation
%
function [fit, ffun, fgrad] = de_ftab_fit_exp(sl, fm, MM, wt, mtype, mac, kev)

sll = ndgrid_jf('mat', sl{:});

if isempty(mac)
	if isempty(mtype), error 'mac or mtype required', end
	mac = xray_read_atten(mtype, kev);
end

LL = length(sl);

if isempty(wt), wt = num2cell(ones(MM,1)); end

Ab = 1;
for ll=1:LL
	sl = col(stackpick(sll, ll));
	Ab = Ab .* exp(-sl * mac(:,ll)'); % [#s*, #E]
end

fit.kev = cell(1,MM);
fit.mac = cell(1,MM);
fit.coef = cell(1,MM);
for mm=1:MM
	if LL == 1
		dat = fm(:,mm);
	elseif LL == 2
		dat = fm(:,:,mm);
	else
		fail 'not done'
	end
	dat = stackpick(fm,mm);
	y = exp(-dat);
	Wh = spdiag(sqrt(wt{mm}(:)), 'nowarn');
	x = wls_simplex(Ab, y(:), Wh); % coefficients for each energy

	ie = x > 1e-6; % find key energies
	fit.kev{mm} = kev(ie);
	fit.mac{mm} = mac(ie,:); % [#E,L]
	A = 1;
	for ll=1:LL
		sl = col(stackpick(sll,ll));
		A = A .* exp(-sl * fit.mac{mm}(:,ll)'); % [#s*, #E]
	end
	fit.coef{mm} = wls_simplex(A, y(:), Wh); % refit with key energies
end
ffun = @de_ftab_fit_exp_eval;
fgrad = @de_ftab_fit_exp_grad;

end % de_ftab_fit_exp()


%
% de_ftab_fit_exp_eval()
% evaluate 
% in
%	sll	[(Nd),L]	stackup of s1,s2,...,s_L
% out
%	f	[(Nd),M]	stackup of f1,f2,...,f_M
%
function f = de_ftab_fit_exp_eval(fit, sll)
Nd = size(sll); LL = Nd(end); Nd = Nd(1:end-1);
sll = reshape(sll, [], LL); % [*Nd,L]
MM = fit.MM;
f = zeros(prod(Nd),MM);
for mm=1:MM
	A = 1;
	mac = fit.mac{mm};
	for ll=1:LL
		sl = sll(:,ll);
		A = A .* exp(-sl * mac(:,ll)'); % [*Nd,ne]
	end
	tmp = -log(A * fit.coef{mm}); % [*Nd,1]
	f(:,mm) = tmp;
end
f = reshape(f, [Nd MM]);

end % de_ftab_fit_exp_eval()


%
% de_ftab_fit_exp_grad()
% evaluate gradient of f for each of the given s vectors.
% in
%	sll	[(Nd),L]	stackup of s1,s2,...,s_L
% out
%	g	[(Nd),L,M]	stackup of gradients of f(s)
%
function g = de_ftab_fit_exp_grad(fit, sll)
Nd = size(sll); LL = Nd(end); Nd = Nd(1:end-1);
sll = reshape(sll, [], LL); % [*Nd,L]
MM = fit.MM;
g = zeros(prod(Nd), LL, MM);
for mm=1:fit.MM
	A = 1;
	mac = fit.mac{mm}; % [ne,L]
	alf = fit.coef{mm}; % [ne,1]
	for ll=1:LL
		sl = sll(:,ll);
		A = A .* exp(-sl * mac(:,ll)'); % [*Nd,ne]
	end
	vm = A * alf; % [*Nd,1]
	tmp = A * (mac .* repmat(alf, [1 LL])); % [*Nd,L]
	g(:,:,mm) = tmp ./ repmat(vm, [1 LL]);
end
g = reshape(g, [Nd LL MM]);

end % de_ftab_fit_exp_grad()


%
% de_ftab_fit_show_sp()
% compare true spectra to fitted spectra
%
function out = de_ftab_fit_show_sp(fit, en, sp)
if nargin < 3, fail 'de_ftab_fit_show_sp(fit, en, sp)', end

if ~streq(fit.type, 'exp3'), printm 'show_sp only done for exp3', return, end
if im
	clf, pl = (fit.MM+1)*100 + 10 + 1;
	subplot(pl)
	plot(en, sp * diag(1 ./ max(sp)))
	if isfield(fit, 'kev')
		for mm=1:fit.MM
			subplot(pl+mm)
			bar(fit.kev{mm}, fit.coef{mm})
			axis tight, axisx(minmax(en))
		end
	else
		warn 'kev unknown'
	end
end

if nargout, out = []; end

end % de_ftab_fit_show_sp()


%
% de_ftab_fit_show_fm()
% compare fit to sampled fm
%
function out = de_ftab_fit_show_fm(fit, sl, fm)
if nargin < 3, error 'de_ftab_fit_show_fm(fit, sl, fm)', end

sll = ndgrid_jf('mat', sl{:});

if fit.MM ~= 2, error 'only M=2 done', end

fh = fit.fmfun(sll);

switch length(sl) % LL
case 1
	s1 = sl{1};
	im clf
	plot(	s1, fm(:,1), 'c.', s1, fm(:,2), 'y.', ...
		s1, fh(:,1), 'c-', s1, fh(:,2), 'y-')
	xlabel 's1'
	ylabel 'fm(s1)'
	legend('true m=1', 'true m=2', 'fit m=1', 'fit m=2', 4)

case 2
	s1 = sl{1};
	s2 = sl{2};
	smax(1) = max(s1);
	smax(2) = max(s2);
	fmax = max(fm(:));

	ax = [0 smax(1) 0 smax(2) 0 fmax];
	cax = [0 fmax];

	im clf
	im pl 2 2

	show(1, s1, s2, fm(:,:,1), ax, cax, 'f_1(s)')
	% text(-60, 10, '[cm^2/g]')
	show(2, s1, s2, fm(:,:,2), ax, cax, 'f_2(s)')

	show(3, s1, s2, fh(:,:,1), ax, cax, 'f_1 approx')
	show(4, s1, s2, fh(:,:,2), ax, cax, 'f_2 approx')
otherwise
	fail 'not done'
end

if nargout, out = []; end

end % de_ftab_fit_show_fm()


%
% de_ftab_fit_show_err()
% show fit errors, and relate to HU
% f = mac s, mac(70kev,H2O) = 0.2 cm^2 / g = 1000 HU
% mh = mac * 1000 HU / (0.2 g/cm^2)
% so Df/Dmh = Df/Dmac * Dmac/Dmh = (50 g/cm^2) (0.2 cm^2/g) / 1000 HU = 1/100 HU
% so Dmh = 100 HU * Df for 50 cm of H2O.
%
function out = de_ftab_fit_show_err(fit, sl, fm)
if nargin < 3, error 'de_ftab_fit_show_err(fit, sl, fm)', end

sll = ndgrid_jf('mat', sl{:});

fh = fit.fmfun(sll);
err = fh - fm;
printm('worst model error = %g of %g', max(abs(err(:))), max(fm(:)))
printm('worst model error = %g HU over 50cm H2O', 100*max(abs(err(:))))
for mm=1:fit.MM
	ee = stackpick(err,mm);
	printm('worst error (mm=%d): %g', mm, max(abs(ee(:))))
end

if max(abs(err(:))) == 0, return, end
if fit.MM > 2, error 'only M=2 done', end
if 0
	e1 = err(:,:,1);
	e2 = err(:,:,2);
	disp([minmax(e1); minmax(e2)]')
	printm('worst error1 %g', max(col(abs(err(:,1,:)))))
	printm('worst error2 %g', max(col(abs(err(:,2,:)))))
end
%err = abs(err);

switch length(sl) % LL
case 1
	s1 = sl{1};
	plot(s1, err(:,1), 'c-', s1, err(:,2), 'y-')
	xlabel 's1'
	ylabel 'error'
	legend('m=1', 'm=2', 4)

case 2
	s1 = sl{1};
	s2 = sl{2};
	smax(1) = max(s1);
	smax(2) = max(s2);

	elim = minmax(err(:))';
	%elim = [-1 1] * 0.01; % +/- 1 HU
	ax = [0 smax(1) 0 smax(2) elim];
	im plc 1 2
	for mm=1:fit.MM
		show(mm, s1, s2, err(:,:,mm), ax, elim, sprintf('f_%d error', mm))
	end
otherwise
	fail 'not done'
end

if nargout, out = []; end

end % de_ftab_fit_show_err()


%
% show()
%
function show(pl, x, y, f, ax, cax, ti)
if ~im, return, end
im('subplot', pl)
if 1
	mesh(x,y,f')
	colormap hsv, caxis(cax), cbar
	axis(ax)
	xtick, ytick, ztick, zwhite
else
	im(x,y,f), cbar
	xtick, ytick
end
xlabel 's_1', ylabel 's_2', title(ti)

end % show()


%
% de_ftab_fit_exp3_test()
%
function de_ftab_fit_exp3_test

%stype = 'mono,70';
stype = 'ps1';
xrs = xray_read_spectra(stype);

if 1 % L=1 component
	sl{1} = linspace(0, 50, 26);
	mtype = {'water'};
	mas = xray_read_mac(mtype);
	mac = mas.mac(xrs.en);
	if im
		clf, semilogy(xrs.en, mac), legend(mtype{:})
	end
	sll = ndgrid_jf('mat', sl{:});
	fm = de_ftab_fm(sll, mac, xrs.Ide);
	fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype);
	g = fit.fgrad(sll);
%	fit = de_ftab_fit(sl, fm, 'type', 'poly')

	if 1
		fit.show_sp(xrs.en, xrs.sp);
		prompt
		fit.show_fm(sl, fm);
		prompt
		fit.show_err(sl, fm);
		prompt
	end
end

if 1 % L=2
	sl{1} = linspace(0, 50, 26);
	sl{2} = linspace(0, 30, 31);
	mtype = {'water', 'bone'};
	mas = xray_read_mac(mtype);
	mac = mas.mac(xrs.en);
	if im
		clf, semilogy(xrs.en, mac), legend(mtype{:})
	end
	sll = ndgrid_jf('mat', sl{:});
	fm = de_ftab_fm(sll, mac, xrs.Ide);
	fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype);
	g = fit.fgrad(sll);
	%fit = de_ftab_fit(sl, fm, 'type', 'poly')

	if 1
		fit.show_sp(xrs.en, xrs.sp);
		prompt
		fit.show_fm(sl, fm);
		prompt
		fit.show_err(sl, fm);
		prompt
	end
end

end % de_ftab_fit_exp3_test()
