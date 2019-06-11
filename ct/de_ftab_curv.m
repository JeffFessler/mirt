  function fit = de_ftab_curv(fit, varargin)
%|function fit = de_ftab_curv(fit, [options])
%|
%| Augment a fitted DE BH function strum object
%| to also have curvature methods for LS (or WLS) data-fit term
%|
%| in
%|	fit	strum		from de_ftab_fit (type: 'exp')
%|
%| option
%|	'ctype'	char		curvature type: (default: 'pre10')
%|				'pre10': precomputed based on 1st deriv at 0
%|				'pre2': precomputed based on surrogates
%|				'newt': newton curvature
%|
%|	'ltype'	char		log-likelihood type: (default: 'ls')
%|				'ls' least-squares sum_m wm (hfim - fm(s))^2 / 2
%|	'show'	0|1		default: 0
%|
%| out
%|	fit	strum		strum object with additional methods:
%|	.ls_cost(fh [1 M], ss [() L], ?wm [1 M])	[()] WLS cost
%|	.ls_grad(fh [1 M], ss [() L], ?wm [1 M])	[() L] WLS gradient
%|	.ls_hess(fh [1 M], ss [() L], ?wm [1 M])	[() L L] WLS hessian
%|	.ls_curv(sll, hfim, fim(sll)|[], ?'ctype', ctype)
%|				 	[() M] curvatures for WLS case
%|	.show_sur()			show surrogates
%|
%| Copyright 2008-10-5, Jeff Fessler, University of Michigan

if nargin == 1 && streq(fit, 'test'), de_ftab_curv_test, clear, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.show = false;
arg.ctype = 'pre10';
arg.ltype = 'ls';
arg = vararg_pair(arg, varargin);

if ~streq(arg.ltype, 'ls'), fail 'only ''ls'' done', end
new.ltype = arg.ltype;
new.ctype = arg.ctype;

switch fit.type
case 'exp'

	if fit.LL == 1
		new.cs = de_ftab_curv_exp_ls_setup(fit, 'show', arg.show);
	end

	g0 = fit.fgrad(zeros(1, fit.LL)); % [1 L M]
	new.cs.g0norm2 = squeeze(sum(g0.^2,2)); % [M]

otherwise
	fail('unknown fit type %s', arg.type)
end

meth = {'show_sur', @de_ftab_curv_show_sur, '()';
	'ls_cost', @de_ftab_curv_exp_ls_cost, ...
		'(fh [1 M], ss [() L], ?wm [1 M]) -> [()]';
	'ls_grad', @de_ftab_curv_exp_ls_grad, ...
		'(fh [1 M], ss [() L], ?wm [1 M]) -> [() L]';
	'ls_hess', @de_ftab_curv_exp_ls_hess, ...
		'(fh [1 M], ss [() L], ?wm [1 M]) -> [() L L]';
	'ls_curv', @de_ftab_curv_exp_ls_curv, ...
		'(sll [() L], hfim [() M], fim(sll)|[] [() M], [''ctype'', ctype]) -> [() M]';
	'ls_surr', @de_ftab_curv_exp_ls_surr, ...
		'(s0 [1 L], fh [1 M], curv [1 M], ss [() L]) -> [() M]';
	};

fit = strum(new, meth, 'base', fit);

if arg.show && im
	fit.show_sur;
end

end % de_ftab_curv()


% de_ftab_curv_exp_ls_setup()
% tabulate maximum of (f'(s))^2 + (-d^/ds^2 f(s)) (hfim - f(s)),
% needed for surrogate curvature for LS cost function
function cs = de_ftab_curv_exp_ls_setup(fit, varargin)

cs.show = 0;
cs.hfim_list = linspace(-5, 10, 31);
cs = vararg_pair(cs, varargin);
cs.show = cs.show & im;

LL = fit.LL;
if LL~=1, fail 'only LL=1 done', end
sl = linspace(-5, 5, 201)';
%sl = ndgrid_jf('mat', sl{:});
f = fit.fmfun(sl); % [*Ns M]
g = fit.fgrad(sl);
g = squeeze(g); % [*Ns M]
h = fit.fhess(sl);
h = squeeze(h); % [*Ns M]

cs.ctable = zeros(length(cs.hfim_list), fit.MM);
warned = 0;
for ii = 1:length(cs.hfim_list)
	hfim = cs.hfim_list(ii);
	tmp = g.^2 - h .* (hfim - f); % [*Ns M]

	[cmax where] = max(tmp); % [1 M]
	if cs.show
		plot(sl, tmp, '-', sl(where), cmax, 'o')
	drawnow
%	prompt
	end

	% if special ntype is used, then cmax should be in middle
	if ~isempty(fit.ntype) && ~warned ...
		&& (any(where == 1) || any(where == length(sl))) && hfim > 0
		warn('cmax not in middle')
		warned = 1;
%		keyboard
	end

	cs.ctable(ii,:) = cmax;
end
if cs.show, prompt, end

% curvature function: independent of sil because it is "max"
cs.curvf = @(cs, mm, sil, hfim) ...
	interp1(cs.hfim_list, cs.ctable(:,mm), hfim, 'linear', 'extrap');

if cs.show
	fh = linspace(-6, 12, 101);
	plot(cs.hfim_list, cs.ctable, 'o', fh, cs.curvf(cs, 1, 0, fh), '-')
	xlabel 'hat f im', ylabel 'curv2'
prompt
end

return

% old way:

tmp = abs(g.^2 + h .* f);
[cmax where] = max(tmp); % [1 M]
if cs.show
	plot(sl, tmp, '-', sl(where), cmax, 'o')
prompt
end
if any(where == 1) || any(where == length(sl))
	fail('cmax not in middle')
end

tmp = abs(h);
[vmax where] = max(tmp); % [1 M]
if any(where == 1) || any(where == length(sl))
	fail('vmax not in middle')
end
if cs.show
	plot(sl, tmp, '-', sl(where), vmax, 'o')
end

curv_max = {cmax, vmax};

end % de_ftab_curv_exp_ls_setup()


% de_ftab_curv_exp_ls_curv()
% surrogate curvatures for LS cost function for each of the given s vectors.
% in
%	fit	strum
%	sil	[(nd) L]	component material density integrals
%	hfim	[(nd) M]	\hat{f}_im (log estimates)
%	fim	[(nd) M]	f_im(sll) - or empty
% option
%	'ctype'
% out
%	curv	[(nd) M]	curvatures for each term of LS cost function
function curv = de_ftab_curv_exp_ls_curv(fit, sil, hfim, fim, varargin)
LL = fit.LL;
MM = fit.MM;
nd = size(hfim); if MM ~= nd(end), error 'hfim size', end; nd = nd(1:end-1);
sil = reshapee(sil, [], LL); % [*nd L]
hfim = reshapee(hfim, [], MM); % [*nd M]
if size(sil,1) ~= prod(nd), error 'sil size', end
curv = zeros(prod(nd), MM); % [*nd M]

arg.ctype = fit.ctype;
arg = vararg_pair(arg, varargin);

switch arg.ctype

case 'newt' % newton curvature
	ls_hess_m = de_ftab_curv_exp_ls_hess_m(fit, hfim, sil); % [L L M *ns]
	if LL == 1
		curv = reshape(ls_hess_m, [MM prod(nd)])'; % [*nd M]
	else
		for is=1:prod(nd)
			for mm=1:MM
				hh = ls_hess_m(:,:,mm,is);
				curv(is,mm) = norm(hh);
			end
		end
	end

case 'pre10'
	curv = fit.cs.g0norm2; % [M 1]
	curv = repmat(curv', [nd 1]); % [(nd) M]

case 'pre2'
%	curv2 = @de_ftab_curv_exp_ls_curv2; % todo: clean up

	%curv_max = fit.curv_max;
	%cmax = curv_max{1}; % [1 M]
	%vmax = curv_max{2}; % [1 M]
	for mm=1:fit.MM
	%	mac = fit.mac{mm}; % [Ne L]
	%	vmax = (max(mac) - min(mac))^2 / 4;
	%	curv(:,mm) = max(mac)^2 + vmax * hfim(:,mm) + cmax(mm);
	%	curv(:,mm) = vmax(mm) * abs(hfim(:,mm)) + cmax(mm);
		curv(:,mm) = fit.cs.curvf(fit.cs, mm, sil, hfim(:,mm));
	end

otherwise
	fail('unknown ctype %s', arg.ctype)
end

curv = reshape(curv, [nd MM]);

end % de_ftab_curv_exp_ls_curv()


% de_ftab_curv_exp_ls_cost()
% evaluate WLS cost function: \sum_m w_m (hf_m - fm(s))^2 / 2
% in
%	fit	strum
%	fh	[1 M]		\hat{f}_im (log estimates)
%	ss	[(ns) L]	component material density integrals
% option
%	wm	[1 M]		weights
% out
%	cost	[(ns) 1]	WLS cost function
%
function cost = de_ftab_curv_exp_ls_cost(fit, fh, ss, wm)
M = fit.MM;
ns = size(ss); if fit.LL ~= ns(end), error 'ss size', end; ns = ns(1:end-1);
ss = reshapee(ss, [], fit.LL); % [*ns L]
fs = fit.fmfun(ss); % [*ns M]
if ~isvar('wm') || isempty(wm)
	wm = ones(1,M);
end
cost = 0;
for mm=1:M
	cost = cost + wm(mm) * (fs(:,mm) - fh(mm)).^2 / 2; % [*ns 1]
end
cost = reshape(cost, [ns 1]); % [(ns)]

end % de_ftab_curv_exp_ls_cost()


% de_ftab_curv_exp_ls_grad()
% gradient of WLS cost function
% in
%	fit	strum
%	fh	[1 M]		\hat{f}_im (log estimates)
%	ss	[(ns) L]	component material density integrals
% option
%	wm	[1 M]		weights
% out
%	grad	[(ns) L]	gradient of WLS cost function
%
function grad = de_ftab_curv_exp_ls_grad(fit, fh, ss, wm)
M = fit.MM;
ns = size(ss); if fit.LL ~= ns(end), error 'ss size', end; ns = ns(1:end-1);
ss = reshapee(ss, [], fit.LL); % [*ns L]
fs = fit.fmfun(ss); % [*ns M]
fgrad = fit.fgrad(ss); % [*ns L M]
if ~isvar('wm') || isempty(wm)
	wm = ones(1,M);
end
grad = 0;
for mm=1:M
	err = fs(:,mm) - fh(mm); % [*ns 1]
	err = repmat(err, [1 fit.LL]); % [*ns M]
	grad = grad + wm(mm) * err .* fgrad(:,:,mm);
end
grad = reshape(grad, [ns fit.LL]); % [(ns) L]

end % de_ftab_curv_exp_ls_grad()


% de_ftab_curv_exp_ls_hess_m()
% hessians of each term of LS cost function
% in
%	fit	strum
%	fh	[1 M]		\hat{f}_im (log estimates)
%	ss	[(ns) L]	component material density integrals
% out
%	hess	[L L M (ns)]	hessians of LS cost function terms
%
function hess = de_ftab_curv_exp_ls_hess_m(fit, fh, ss)
L = fit.LL;
M = fit.MM;
ns = size(ss); if L ~= ns(end), error 'ss size', end; ns = ns(1:end-1);
ss = reshapee(ss, [], L); % [*ns L]
fs = fit.fmfun(ss); % [*ns M]
fgrad = fit.fgrad(ss); % [*ns L M]
fhess = fit.fhess(ss); % [*ns L L M]
hess = zeros([L L M prod(ns)]);
for mm=1:M
	err = fs(:,mm) - fh(mm); % [*ns 1]
	if L == 1
		h1 = fgrad(:,1,mm).^2; % [*ns 1]
		h2 = err .* fhess(:,1,1,mm);
		hess(1,1,mm,:) = h1 + h2;
	else
%		err = repmat(err, [1 L]); % [*ns M]
		for is=1:prod(ns)
			gg = squeeze(fgrad(is,:,mm)); % [1 L]
			h1 = gg' * gg; % H1
			tmp = reshape(fhess(is,:,:,mm), L, L);
			h2 = err(is) * tmp;
			hess(:,:,mm,is) = h1 + h2;
		end
	end
end
hess = reshape(hess, [L L M ns]); % [L L M (ns)]
%pr cond(hess)

end % de_ftab_curv_exp_ls_hess_m()


% de_ftab_curv_exp_ls_hess()
% hessians of LS cost function
% in
%	fit	strum
%	fh	[1 M]		\hat{f}_im (log estimates)
%	ss	[(ns) L]	component material density integrals
% option
%	wm	[1 M]		weights
% out
%	hess	[L L (ns)]	hessians of LS cost function
%
function hess = de_ftab_curv_exp_ls_hess(fit, fh, ss, wm)
L = fit.LL;
M = fit.MM;
hess = de_ftab_curv_exp_ls_hess_m(fit, fh, ss); % [L L M (ns)]
if isvar('wm') && ~isempty(wm)
	dim = size(hess);
	dim = dim(4:end);
	hess = reshape(hess, [L L M prod(dim)]); % [L L M *ns]
	tmp = 0;
	for mm=1:M
		tmp = tmp + wm(mm) * hess(:,:,mm,:); % [L L *ns]
	end
	hess = reshape(tmp, [L L dim]); % [L L (ns)]

else
	hess = sum(hess, 3); % [L L 1 (ns)]
	tmp = size(hess);
	hess = reshape(hess, [L L tmp(4:end)]); % [L L (ns)]
end

end % de_ftab_curv_exp_ls_hess()


% de_ftab_curv_exp_ls_surr()
% evaluate surrogate for WLS cost function: \sum_m w_m (hf_m - fm(s))^2 / 2
% in
%	fit	strum
%	s0	[1 L]		expansion point
%	fh	[1 M]		\hat{f}_im (log estimates)
%	curv	[1 M]		curvatures
%	ss	[(ns) L]	component material density integrals
% option
%	wm	[1 M]		weights
% out
%	qs	[(ns) 1]	quadratic surrogate for LS cost function
%
function qs = de_ftab_curv_exp_ls_surr(fit, s0, fh, curv, ss, wm)
ns = size(ss); if fit.LL ~= ns(end), error 'ss size', end; ns = ns(1:end-1);
ss = reshapee(ss, [], fit.LL); % [*ns L]

if ~isvar('wm') || isempty(wm)
	wm = ones(1, fit.MM);
end

f0 = fit.fmfun(s0); % [1 M]

g0 = fit.fgrad(s0); % [1 L M]
g0 = reshape(g0, fit.LL, fit.MM); % [L M]

% 0th-order term of surrogate
q0 = sum((f0 - fh).^2) / 2; % [1]

% 1st-order term of surrogate
sdif = ss - repmat(s0, [prod(ns) 1]); % [*ns L]
q1 = sdif * g0 * diag(wm) * (f0 - fh)'; % [*ns 1]

% 2nd-order term of surrogate
q2 = 1/2 * (curv * wm') * sum(sdif.^2,2); % [*ns 1]

qs = q0 + q1 + q2; % [*ns 1]
qs = reshape(qs, [ns 1]); % [(ns)]

end % de_ftab_curv_exp_ls_surr()


% de_ftab_curv_show_sur1()
% show surrogates for the case L=1
function de_ftab_curv_show_sur1(fit, varargin)

s = linspace(-4,16,201)';
f = fit.fmfun(s); % [ns M]

hf1list = repmat([0 1 2 4]', [1 fit.MM]);
s1list = [0 1 2 3]' / 0.25;

if im, jf plc 2 2, end
for ii=1:length(s1list)
	hf1 = hf1list(ii,:);
	s1 = s1list(ii);
	f1 = fit.fmfun(s1);
	g1 = fit.fgrad(s1);
	g1 = squeeze(g1)';

	clist = {'pre2', 'pre10', 'newt'};

	for mm=1:fit.MM
		cost = (hf1(:,mm)-f(:,mm)).^2 / 2;
		for ic=1:length(clist)
			curv = fit.ls_curv(s1, hf1, f1, 'ctype', clist{ic});
			q1(:,ic) = (hf1(mm) - f1(mm))^2/2 ...
				+ g1(mm) * (f1(mm)-hf1(mm)) * (s-s1) ...
				+ curv(mm)/2 * (s-s1).^2;
		end

		if im
			figure(mm)
			if ii==1, clf, end
			jf('sub', ii)
			plot(s, cost, '-', s, q1, '--')
			axis([minmax(s)' minmax(cost(s >= -1))'])
			if ii==1
				legend('cost', clist{:})
				titlef('mm=%d', mm)
			end
		end
	end
end

prompt
if im, figure(2), close, end

end % de_ftab_curv_show_sur1()


% de_ftab_curv_show_sur2()
% show surrogates for the case L=2
function de_ftab_curv_show_sur2(fit, varargin)

sls = de_ftab_sls('min', [-1 -1], 'max', [20 10], 'n', [21 11]);
%sls = de_ftab_sls('min', [-1 -1], 'max', [2 2], 'n', [201 301]);

s0 = [10 5];
%s0 = [0 0];

f0 = fit.fmfun(s0);
fh = f0 + 0*[0.8 -0.6]; % "noise"

%pr fit.ctype

clist = {'newt', 'pre10'};%, 'pre2'}; % todo later
for ic=1:length(clist)
	ctype = clist{ic};
	curv = fit.ls_curv(s0, fh, f0, 'ctype', ctype); % [1 M]
	pr sum(curv)

	qs = fit.ls_surr(s0, fh, curv, sls.sll); % [(ns)]
	cost = fit.ls_cost(fh, sls.sll);

	if 0 && im
%		qs(qs > max(cost(:))) = nan;
		clim = minmax(cost);
		clf, mesh(sls.sll(:,:,1), sls.sll(:,:,2), cost)
		ax = axis;
		hold on
		mesh(sls.sll(:,:,1), sls.sll(:,:,2), qs)
		hold off
	%	axis(ax)
		caxis(clim)
		colormap hsv, cbar
		xlabel 's1', ylabel 's2'
	return
	end

%	clf, plot3(sls.sll(:,:,1), sls.sll(:,:,2), cost)
%	clf, im(s1, s2, cost)

	h0 = fit.ls_hess(fh, s0); % local hessian
	pr cond(h0)
	[v d] = eig(h0);
	v0 = v(:,imax(abs(diag(d))))'; % highest curvature direction
	drawvect = @(p0, p1) plot([p0(1) p1(1)], [p0(2) p1(2)], '-o');

	if im % contours
		clf
		fun = @(c) sqrt(c);
		v = sls.sll(:,:,1) >= 0 & sls.sll(:,:,2) >= 0;
		v = linspace(0, max(fun(cost(v))), 20);
		contour(sls.sll(:,:,1), sls.sll(:,:,2), fun(cost), v, 'c-')
		axis equal
		hold on
		contour(sls.sll(:,:,1), sls.sll(:,:,2), fun(qs), v, 'y-')
		hold off
		hold on
		drawvect(s0-2*v0, s0+2*v0)
		hold off
	prompt
	end

	if im % plot profile along worst direction
%		t = linspace(-3, 3, 101)';
		t = linspace(min(sls.sll(:)), max(sls.sll(:)), 201)';
%		[v d] = eig(fit.mac_eff)
		ss = (t - mean(s0)) * v0; % [ns L]
		ss = ss + repmat(s0, [length(t) 1]); % [ns L]
		q1 = fit.ls_surr(s0, fh, curv, ss);
		k1 = fit.ls_cost(fh, ss);
		plot(t, k1, 'c-', t, q1, 'y-')
		axisy(minmax(k1(min(ss,[],2) > 0)))
		title('profile along worst directon')
	prompt
	end
		
	if 0 && im, one profile
		s1 = sls.sl{1};
		s2 = sls.sl{2};
		i2 = imin(abs(s2 - s0(2)));
		jf_equal(s2(i2), s0(2))
		plot(s1, cost(:,i2), 'c-x')
		axis tight
		hold on
		plot(s1, qs(:,i2), 'y-o')
		hold off
	end

end % for ic

end % de_ftab_curv_show_sur2()


% de_ftab_curv_show_sur()
function out = de_ftab_curv_show_sur(fit, varargin)

switch fit.LL
case 1
	de_ftab_curv_show_sur1(fit, varargin)
case 2
	de_ftab_curv_show_sur2(fit, varargin)
otherwise
	warn('show L=%d not done', fit.LL)
end

if nargout, out = []; end
end % de_ftab_curv_show_sur()


% de_ftab_curv_test1()
% L=1
function de_ftab_curv_test1

stype = 'ps1';
stype = 'poly1,80,140'; % todo: filter!
xrs = xray_read_spectra(stype);
mtype = 'water';
mas = xray_read_mac(mtype);
sl = {linspace(0, 50, 26)'};

mac = mas.mac(xrs.en);
sll = ndgrid_jf('mat', sl{:});
fm = de_ftab_fm(sll, mac, xrs.Ide);

fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype, 'ntype', '');
fit = de_ftab_curv(fit, 'show', 0);
fit.show_sur;

end % de_ftab_curv_test1()



% de_ftab_curv_test2()
% L=2
function de_ftab_curv_test2

stype = 'poly1,80,140';
%stype = 'mono,60,100'
xrs = xray_read_spectra(stype);
mtype = {'water', 'bone'};
mas = xray_read_mac(mtype);
sls = de_ftab_sls;
fm = de_ftab_fm(sls.sll, mas.mac(xrs.en), xrs.Ide);
fit = de_ftab_fit(sls.sl, fm, 'type', 'exp', 'mtype', mtype, 'ntype', '');
ctype = 'pre10';
ctype = 'newt';
fit = de_ftab_curv(fit, 'ctype', ctype, 'show', 0);
fit.show_sur;

end % de_ftab_curv_test2()


% de_ftab_curv_test()
function de_ftab_curv_test
de_ftab_curv_test2
de_ftab_curv_test1

end % de_ftab_curv_test()
